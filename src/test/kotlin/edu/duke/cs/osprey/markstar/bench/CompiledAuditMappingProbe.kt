package edu.duke.cs.osprey.markstar.bench

import edu.duke.cs.osprey.confspace.SimpleConfSpace
import edu.duke.cs.osprey.confspace.VoxelShape
import edu.duke.cs.osprey.confspace.compiled.ConfSpace as CompiledConfSpace
import edu.duke.cs.osprey.confspace.compiled.PosInterDist
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator
import edu.duke.cs.osprey.energy.compiled.CudaConfEnergyCalculator
import edu.duke.cs.osprey.gpu.Structs
import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.compiler.ConfSpaceCompiler
import edu.duke.cs.osprey.gui.forcefield.Forcefield
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.io.toMolecule
import edu.duke.cs.osprey.gui.io.toTomlKey
import edu.duke.cs.osprey.gui.io.toBytes
import edu.duke.cs.osprey.gui.io.withService
import edu.duke.cs.osprey.gui.motions.DihedralAngle
import edu.duke.cs.osprey.gui.prep.Anchor
import edu.duke.cs.osprey.gui.prep.ConfSpace as PrepConfSpace
import edu.duke.cs.osprey.gui.prep.DesignPosition
import edu.duke.cs.osprey.gui.prep.Proteins
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Polymer
import java.io.File
import java.util.Locale
import org.joml.Vector3d


/**
 * Probe for the long-term compiled-audit path.
 *
 * The legacy audit CSV uses SimpleConfSpace RC indices. The compiled CUDA
 * calculator uses compiled Conf indices, so this tool refuses to treat indices
 * as interchangeable and reports the mapping explicitly.
 */
object CompiledAuditMappingProbe {

	@JvmStatic
	fun main(args: Array<String>) {

		val pdbPath = requireProperty("osprey.bench.pdbPath")
		val mutable = System.getProperty("osprey.bench.mutable", "")
		val flexible = System.getProperty("osprey.bench.flexible", "")
		val designId = System.getProperty("osprey.bench.designId", "unknown")
			val state = System.getProperty("osprey.audit.state", "Complex")
			val doCompile = java.lang.Boolean.getBoolean("osprey.compiledAudit.compile")
			val doCudaSweep = java.lang.Boolean.getBoolean("osprey.compiledAudit.cudaSweep")
			val previewConfs = Integer.getInteger("osprey.compiledAudit.previewConfs", 8)

			println("==============================================")
		println("  Compiled Audit Mapping Probe")
		println("  design=$designId state=$state")
		println("  pdb=$pdbPath")
			println("  mutable=$mutable")
			println("  flexible=$flexible")
			println("  compile=$doCompile")
			println("  cudaSweep=$doCudaSweep")
			println("==============================================")

		val confSpaces = AuditLeafCCD.buildConfSpaces(pdbPath, mutable, flexible)
		val simple = AuditLeafCCD.selectConfSpace(confSpaces, state)
		printLegacy(simple, previewConfs)

		val conflib = ConfLib.from(OspreyGui.getResourceAsString("conflib/lovell.conflib"))
		val prepBuild = buildPrepConfSpace(simple, conflib)
		printPrep(prepBuild, previewConfs)

		if (!doCompile) {
			println("compiled: skipped, set -Dosprey.compiledAudit.compile=true to compile")
			return
		}

		val compiled = try {
			compile(prepBuild.confSpace)
		} catch (e: RuntimeException) {
			println("compiled: failed ${e.message}")
			rootCause(e)?.let { cause ->
				if (cause !== e) {
					println("compiled.rootCause ${cause.javaClass.name}: ${cause.message}")
				}
			}
			return
		}
		printCompiled(compiled, previewConfs)

			val mappings = buildMappings(simple, compiled)
			printMappings(mappings)
			printSampleAssignment(simple, mappings)
			if (doCudaSweep) {
				runCudaSweep(simple, compiled, mappings)
			}
		}

	private fun requireProperty(name: String) =
		System.getProperty(name)?.takeIf { it.isNotBlank() }
			?: throw IllegalArgumentException("missing required system property: $name")

	private data class PrepBuild(
		val confSpace: PrepConfSpace,
		val links: List<PositionLink>,
		val issues: List<String>
	)

	private data class PositionLink(
		val simple: SimpleConfSpace.Position,
		val prep: DesignPosition
	)

	private data class PrepMol(
		val mol: Polymer,
		val strandIndex: Int,
		val fragmentIndex: Int
	)

	private fun buildPrepConfSpace(simple: SimpleConfSpace, conflib: ConfLib): PrepBuild {

		val prepMols = simple.strands.flatMapIndexed { index, strand ->
			val mol = normalizeAmberAtomNames(strand.mol.toMolecule("strand%02d".format(index)) as? Polymer
				?: throw IllegalArgumentException("compiled audit probe expects legacy strands to convert to Polymer")
			)
			splitPeptideConnectedFragments(mol, index)
		}

		val confSpace = PrepConfSpace(prepMols.map { MoleculeType.Protein to it.mol }).apply {
			name = "Compiled audit probe"
			conflibs.add(conflib)
		}

		val links = ArrayList<PositionLink>()
		val issues = ArrayList<String>()
		for (simplePos in simple.positions) {

			val resId = AuditResidueId.parse(simplePos.resNum)
			val prepMol = findPrepMol(prepMols, resId)
			val residue = prepMol.mol.findChainOrThrow(resId.chainId).findResidueOrThrow(resId.resId)
			val prepPos = Proteins.makeDesignPosition(prepMol.mol, residue, simplePos.resNum)
			confSpace.addPosition(prepPos)

			val posConfSpace = confSpace.positionConfSpaces.getOrMake(prepPos)
			for (resType in simplePos.resFlex.resTypes) {
				val type = normalizeType(resType)
				try {
					addLegacyLibraryConformations(conflib, simplePos, prepPos, posConfSpace, resType)
					posConfSpace.mutations.add(type)
				} catch (e: RuntimeException) {
					issues.add("${simplePos.resNum}: cannot add library type $type: ${e.message}; ${describePosition(prepPos)}")
				}
			}
			if (simplePos.resFlex.addWildTypeRotamers) {
				val type = normalizeType(simplePos.resFlex.wildType)
				try {
					val fallback = addLegacyWildTypeConformation(conflib, simplePos, prepPos, posConfSpace)
					posConfSpace.mutations.add(type)
					if (fallback != null) {
						issues.add("${simplePos.resNum}: $fallback")
					}
				} catch (e: RuntimeException) {
					issues.add("${simplePos.resNum}: cannot add wild-type conformation $type: ${e.message}; ${describePosition(prepPos)}")
				}
			}
			try {
				addContinuousMotions(prepPos, posConfSpace)
			} catch (e: RuntimeException) {
				issues.add("${simplePos.resNum}: cannot add continuous motions: ${e.message}")
			}

			links.add(PositionLink(simplePos, prepPos))
		}

		return PrepBuild(confSpace, links, issues)
	}

	private fun normalizeAmberAtomNames(src: Polymer): Polymer {
		val dst = Polymer(src.name)
		dst.type = src.type
		dst.netCharge = src.netCharge

		val atomMap = java.util.IdentityHashMap<Atom, Atom>()
		for (chain in src.chains) {
			val dstChain = Polymer.Chain(chain.id)
			for (res in chain.residues) {
				val dstAtoms = res.atoms.map { atom ->
					val dstAtom = Atom(atom.element, amberAtomName(res.type, atom.name), Vector3d(atom.pos))
					atomMap[atom] = dstAtom
					dst.atoms.add(dstAtom)
					dstAtom
				}
				dstChain.residues.add(Polymer.Residue(res.id, res.type, dstAtoms))
			}
			dst.chains.add(dstChain)
		}

		for (atom in src.atoms) {
			val dstAtom = atomMap[atom] ?: continue
			for (bonded in src.bonds.bondedAtoms(atom)) {
				val dstBonded = atomMap[bonded] ?: continue
				dst.bonds.add(dstAtom, dstBonded)
			}
		}

		return dst
	}

	private fun splitPeptideConnectedFragments(src: Polymer, strandIndex: Int): List<PrepMol> {
		val out = ArrayList<PrepMol>()
		var fragmentIndex = 0

		for (chain in src.chains) {
			for (residues in peptideConnectedComponents(src, chain)) {
				val mol = copyChainFragment(
					src,
					chain.id,
					residues,
					"strand%02d_%s_frag%02d".format(strandIndex, chain.id, fragmentIndex)
				)
				out.add(PrepMol(mol, strandIndex, fragmentIndex))
				fragmentIndex += 1
			}
		}

		if (out.isEmpty()) {
			throw IllegalArgumentException("legacy strand $strandIndex converted to no peptide-connected fragments")
		}
		return out
	}

	private fun peptideConnectedComponents(mol: Polymer, chain: Polymer.Chain): List<List<Polymer.Residue>> {
		val residues = chain.residues
		if (residues.isEmpty()) return emptyList()

		val atomToResidueIndex = java.util.IdentityHashMap<Atom, Int>()
		for ((resi, res) in residues.withIndex()) {
			for (atom in res.atoms) {
				atomToResidueIndex[atom] = resi
			}
		}

		val adjacency = MutableList(residues.size) { java.util.LinkedHashSet<Int>() }
		for ((resi, res) in residues.withIndex()) {
			for (atom in res.atoms) {
				if (!isPeptideBackboneEndpoint(atom)) continue
				for (bonded in mol.bonds.bondedAtoms(atom)) {
					val bondedResi = atomToResidueIndex[bonded] ?: continue
					if (bondedResi == resi) continue
					if (!isPeptideBondAtomPair(atom, bonded)) continue
					adjacency[resi].add(bondedResi)
					adjacency[bondedResi].add(resi)
				}
			}
		}

		val seen = BooleanArray(residues.size)
		val components = ArrayList<List<Polymer.Residue>>()
		for (start in residues.indices) {
			if (seen[start]) continue

			val indices = ArrayList<Int>()
			val queue = java.util.ArrayDeque<Int>()
			seen[start] = true
			queue.add(start)

			while (!queue.isEmpty()) {
				val resi = queue.removeFirst()
				indices.add(resi)
				for (next in adjacency[resi]) {
					if (seen[next]) continue
					seen[next] = true
					queue.add(next)
				}
			}

			indices.sort()
			components.add(indices.map { residues[it] })
		}
		return components
	}

	private fun isPeptideBackboneEndpoint(atom: Atom) =
		atom.name.equals("C", ignoreCase = true) || atom.name.equals("N", ignoreCase = true)

	private fun isPeptideBondAtomPair(a: Atom, b: Atom) =
		(a.name.equals("C", ignoreCase = true) && b.name.equals("N", ignoreCase = true))
			|| (a.name.equals("N", ignoreCase = true) && b.name.equals("C", ignoreCase = true))

	private fun copyChainFragment(
		src: Polymer,
		chainId: String,
		residues: List<Polymer.Residue>,
		name: String
	): Polymer {
		val dst = Polymer(name)
		dst.type = src.type
		dst.netCharge = src.netCharge

		val atomMap = java.util.IdentityHashMap<Atom, Atom>()
		val dstChain = Polymer.Chain(chainId)
		val isSingleResidueFragment = residues.size == 1
		for ((resi, res) in residues.withIndex()) {
			val isFragmentNTerminus = resi == 0
			val dstAtoms = res.atoms.mapNotNull { atom ->
				val atomName = amberFragmentAtomName(res, atom, isFragmentNTerminus, isSingleResidueFragment)
					?: return@mapNotNull null
				val dstAtom = Atom(atom.element, atomName, Vector3d(atom.pos))
				atomMap[atom] = dstAtom
				dst.atoms.add(dstAtom)
				dstAtom
			}
			dstChain.residues.add(Polymer.Residue(res.id, res.type, dstAtoms))
		}
		dst.chains.add(dstChain)

		for (atom in src.atoms) {
			val dstAtom = atomMap[atom] ?: continue
			for (bonded in src.bonds.bondedAtoms(atom)) {
				val dstBonded = atomMap[bonded] ?: continue
				dst.bonds.add(dstAtom, dstBonded)
			}
		}

		return dst
	}

		private fun amberFragmentAtomName(
			res: Polymer.Residue,
			atom: Atom,
			isFragmentNTerminus: Boolean,
			isSingleResidueFragment: Boolean
		): String? {
			if (isSingleResidueFragment) {
				if (
					atom.name.equals("H1", ignoreCase = true)
					&& res.atoms.none { it !== atom && it.name.equals("H", ignoreCase = true) }
				) {
					return "H"
				}
				if (atom.name.equals("H2", ignoreCase = true) || atom.name.equals("H3", ignoreCase = true)) {
					return null
				}
			}
			if (
				!isSingleResidueFragment
				&&
				isFragmentNTerminus
				&& atom.name.equals("H", ignoreCase = true)
				&& res.atoms.none { it !== atom && it.name.equals("H1", ignoreCase = true) }
			) {
			return "H1"
		}
		return atom.name
	}

	private fun amberAtomName(resType: String, atomName: String) =
		when {
			// Legacy OSPREY C-terminal CYS templates can name the thiol hydrogen HSG,
			// but Amber/LEaP emits and maps the standard CYS sulfur hydrogen as HG.
			resType.equals("CYS", ignoreCase = true) && atomName.equals("HSG", ignoreCase = true) -> "HG"
			else -> atomName
		}

	private fun addLegacyLibraryConformations(
		conflib: ConfLib,
		simplePos: SimpleConfSpace.Position,
		prepPos: DesignPosition,
		posConfSpace: PrepConfSpace.PositionConfSpace,
		resType: String
	) {
		val legacyTypes = legacyLibraryRawTypes(simplePos, resType)
		if (legacyTypes.isEmpty()) {
			throw IllegalArgumentException("legacy SimpleConfSpace has no library RCs for requested type ${normalizeType(resType)}")
		}

		var added = 0
		val failures = ArrayList<String>()
		for (legacyType in legacyTypes) {
			val candidates = libraryFragmentsForLegacyType(conflib, legacyType, resType)
			if (candidates.isEmpty()) {
				failures.add("$legacyType: no Lovell fragments")
				continue
			}
			val compatible = candidates.filter { prepPos.isFragmentCompatible(it) }
			if (compatible.isEmpty()) {
				failures.add("$legacyType: no compatible fragments among ${candidates.map { it.id }}")
				continue
			}
			for (frag in compatible) {
				posConfSpace.confs.addAll(frag)
				added += frag.confs.size
			}
		}

		if (added <= 0) {
			throw IllegalArgumentException(failures.joinToString("; "))
		}
	}

	private fun addLegacyWildTypeConformation(
		conflib: ConfLib,
		simplePos: SimpleConfSpace.Position,
		prepPos: DesignPosition,
		posConfSpace: PrepConfSpace.PositionConfSpace
	): String? {

		val wtType = legacyWildTypeRawType(simplePos)
		val motionFrags = libraryFragmentsForLegacyType(conflib, wtType, wtType)
			.filter { prepPos.isFragmentCompatible(it) }

		val motionFailures = ArrayList<String>()
		for (motionFrag in motionFrags) {
			try {
				addWildTypeFragment(prepPos, posConfSpace, motionFrag.motions)
				return null
			} catch (e: RuntimeException) {
				motionFailures.add("${motionFrag.id}: ${e.message}")
			}
		}

		// Diagnostic fallback: this restores the WT RC count for mapping work, but is not yet
		// energy-equivalent because the WT RC has no copied library motions.
		addWildTypeFragment(prepPos, posConfSpace, emptyList())
		return "wild-type conformation ${normalizeType(wtType)} added without copied library motions; motion copy failures=${motionFailures.ifEmpty { listOf("no compatible motion fragment") }.joinToString("; ")}"
	}

	private fun addWildTypeFragment(
		prepPos: DesignPosition,
		posConfSpace: PrepConfSpace.PositionConfSpace,
		motions: List<ConfLib.ContinuousMotion>
	) {
		val wtFrag = prepPos.makeFragment(
			"wt-${prepPos.mol.name.toTomlKey()}-${prepPos.name.toTomlKey()}",
			"WildType @ ${prepPos.mol.name} ${prepPos.name}",
			"conf1",
			"conf1",
			motions
		)
		posConfSpace.wildTypeFragment = wtFrag
		posConfSpace.confs.addAll(wtFrag)
	}

	private fun legacyLibraryRawTypes(simplePos: SimpleConfSpace.Position, resType: String): List<String> {
		val normalized = normalizeType(resType)
		return simplePos.resConfs
			.filter { it.type == SimpleConfSpace.ResidueConf.Type.Library }
			.map { it.template.name.toUpperCase(Locale.US) }
			.distinct()
			.filter { normalizeType(it) == normalized }
	}

	private fun legacyWildTypeRawType(simplePos: SimpleConfSpace.Position): String {
		return simplePos.resConfs
			.firstOrNull { it.type == SimpleConfSpace.ResidueConf.Type.WildType }
			?.template
			?.name
			?.toUpperCase(Locale.US)
			?: simplePos.resFlex.wildType.toUpperCase(Locale.US)
	}

	private fun libraryFragmentsForLegacyType(conflib: ConfLib, legacyType: String, requestedType: String): List<ConfLib.Fragment> {
		val legacy = when (legacyType.toUpperCase(Locale.US)) {
			// Legacy Amber configs treat HIS as the doubly-protonated histidine; GUI Lovell
			// stores that concrete protonation variant as HIP, with fragment type HIS.
			"HIS" -> "HIP"
			// Amber/pdb4amber uses CYX for disulfide-bonded cysteine. Lovell has CYS
			// rotamers, and the existing disulfide bond remains in the prep molecule.
			"CYX" -> "CYS"
			else -> legacyType.toUpperCase(Locale.US)
		}
		val exact = conflib.fragments[legacy]
		if (exact != null && normalizeType(exact.type) == normalizeType(requestedType)) {
			return listOfNotNull(
				exact,
				conflib.fragments["${legacy}n"]?.takeIf { normalizeType(it.type) == normalizeType(requestedType) }
			)
		}
		return conflib.fragments.values
			.filter { normalizeType(it.type) == normalizeType(requestedType) }
			.sortedBy { it.id }
	}

	private fun findPrepMol(prepMols: List<PrepMol>, resId: AuditResidueId): PrepMol {
		return prepMols.firstOrNull { prepMol ->
			prepMol.mol.chains.any { chain ->
				chain.id == resId.chainId && chain.residues.any { it.id == resId.resId }
			}
		} ?: throw NoSuchElementException("no legacy strand contains residue ${resId.chainId}${resId.resId}")
	}

	private data class AuditResidueId(val chainId: String, val resId: String) {
		companion object {
			fun parse(resNum: String): AuditResidueId {
				if (resNum.length < 2) {
					throw IllegalArgumentException("expected chain-prefixed residue id, got: $resNum")
				}
				return AuditResidueId(resNum.substring(0, 1), resNum.substring(1).trim())
			}
		}
	}

	private fun normalizeType(type: String) =
		when (type.toUpperCase(Locale.US)) {
			"HID", "HIE", "HIP" -> "HIS"
			"CYX" -> "CYS"
			else -> type.toUpperCase(Locale.US)
		}

	private fun addContinuousMotions(pos: DesignPosition, posConfSpace: PrepConfSpace.PositionConfSpace) {
		val settings = DihedralAngle.LibrarySettings(
			radiusDegrees = VoxelShape.DefaultHalfWidthDegrees,
			includeHydroxyls = true,
			includeNonHydroxylHGroups = true
		)
		for (conf in posConfSpace.confs) {
			conf.motions.addAll(DihedralAngle.ConfDescription.makeFromLibrary(pos, conf.frag, conf.conf, settings))
		}
	}

	private fun compile(confSpace: PrepConfSpace): CompiledConfSpace = withService {
		val compileThreads = Math.max(1, Integer.getInteger("osprey.compiledAudit.compileThreads", 1))
		println("compiled: compileThreads=$compileThreads")
		val bytes = ConfSpaceCompiler(confSpace).run {
			forcefields.add(Forcefield.Amber96)
			forcefields.add(Forcefield.EEF1.configure {
				scale = 1.0
			})
			compile(compileThreads).run {
				waitForFinish()
				report?.compiled?.toBytes()
					?: report?.error?.let { throw IllegalStateException("compiled confspace failed: $it") }
					?: throw IllegalStateException("compiled confspace failed without a report")
			}
		}
		CompiledConfSpace.fromBytes(bytes)
	}

	private fun rootCause(t: Throwable): Throwable? {
		var cause: Throwable = t
		while (cause.cause != null && cause.cause !== cause) {
			cause = cause.cause!!
		}
		return cause
	}

	private fun printLegacy(simple: SimpleConfSpace, previewConfs: Int) {
		println("legacy: positions=${simple.positions.size} confSpaceSize=${simple.numConformations}")
		for (pos in simple.positions) {
			println("legacy.pos index=${pos.index} name=${pos.resNum} mutable=${pos.hasMutations()} rcCount=${pos.resConfs.size}")
			println("  legacy.types ${pos.resConfs.groupBy { normalizeType(it.template.name) }.mapValues { it.value.size }}")
			println("  legacy.rawTypes ${pos.resConfs.groupBy { it.template.name }.mapValues { it.value.size }}")
			for (rc in pos.resConfs.take(previewConfs)) {
				println("  legacy.rc index=${rc.index} type=${normalizeType(rc.template.name)} rawType=${rc.template.name} kind=${rc.type} rot=${rc.rotamerIndex} code=${rc.rotamerCode} dofs=${rc.dofBounds.size}")
			}
		}
	}

	private fun printPrep(build: PrepBuild, previewConfs: Int) {
		val confSpace = build.confSpace
		println("prep: molecules=${confSpace.mols.size} positions=${confSpace.positions().size} confSpaceSize=${confSpace.positionConfSpaces.confSpaceSize()}")
		for ((type, mol) in confSpace.mols) {
			val chains = (mol as? Polymer)
				?.chains
				?.joinToString(";") { chain ->
					val first = chain.residues.firstOrNull()?.id ?: "-"
					val last = chain.residues.lastOrNull()?.id ?: "-"
					"${chain.id}:${chain.residues.size}[$first-$last]"
				}
				?: "nonpolymer"
			println("prep.mol name=${mol.name} type=$type atoms=${mol.atoms.size} chains=$chains")
		}
		for (issue in build.issues) {
			println("prep.issue $issue")
		}
		for (link in build.links) {
			val posConfSpace = confSpace.positionConfSpaces.getOrMake(link.prep)
			val confs = posConfSpace.confs.toList()
			println("prep.pos legacy=${link.simple.resNum} name=${link.prep.name} mutations=${posConfSpace.mutations.sorted()} confCount=${confs.size}")
			println("  prep.position ${describePosition(link.prep)}")
			println("  prep.types ${confs.groupBy { normalizeType(it.frag.type) }.mapValues { it.value.size }}")
			println("  prep.frags ${confs.groupBy { it.frag.id }.mapValues { it.value.size }}")
			for ((index, conf) in confs.take(previewConfs).withIndex()) {
				println("  prep.conf ordinal=$index id=${conf.frag.id}:${conf.conf.id} type=${normalizeType(conf.frag.type)} motions=${conf.motions.size}")
			}
		}
	}

	private fun printCompiled(compiled: CompiledConfSpace, previewConfs: Int) {
		println("compiled: name=${compiled.name} positions=${compiled.positions.size} forcefields=${compiled.forcefieldIds.joinToString("+")}")
		for (pos in compiled.positions) {
			println("compiled.pos index=${pos.index} name=${pos.name} wt=${normalizeType(pos.wildType)} confCount=${pos.confs.size}")
			println("  compiled.types ${pos.confs.groupBy { normalizeType(it.type) }.mapValues { it.value.size }}")
			for (conf in pos.confs.take(previewConfs)) {
				println("  compiled.conf index=${conf.index} id=${conf.id} type=${normalizeType(conf.type)} motions=${conf.motions.size}")
			}
		}
	}

	private data class PositionMapping(
		val simplePos: SimpleConfSpace.Position,
		val compiledPos: CompiledConfSpace.Pos?,
		val rcToConf: Map<Int, Int>,
		val issues: List<String>
	) {
		val complete get() = compiledPos != null && issues.isEmpty() && rcToConf.size == simplePos.resConfs.size
	}

	private fun buildMappings(simple: SimpleConfSpace, compiled: CompiledConfSpace): List<PositionMapping> {
		return simple.positions.map { simplePos ->
			val compiledPos = compiled.positions.find { it.name == simplePos.resNum }
			if (compiledPos == null) {
				return@map PositionMapping(simplePos, null, emptyMap(), listOf("missing compiled position named ${simplePos.resNum}"))
			}

			val rcToConf = LinkedHashMap<Int, Int>()
			val issues = ArrayList<String>()

			val simpleLibByType = simplePos.resConfs
				.filter { it.type == SimpleConfSpace.ResidueConf.Type.Library }
				.groupBy { normalizeType(it.template.name) }
			val compiledLibByType = compiledPos.confs
				.filter { !it.id.startsWith("wt-") }
				.groupBy { normalizeType(it.type) }

			for ((type, rcs) in simpleLibByType) {
				val confs = compiledLibByType[type] ?: emptyList()
				if (rcs.size == confs.size) {
					for ((rc, conf) in rcs.sortedBy { it.index }.zip(confs.sortedBy { it.index })) {
						rcToConf[rc.index] = conf.index
					}
				} else {
					issues.add("library type $type count mismatch legacy=${rcs.size} compiled=${confs.size}")
				}
			}

			val simpleWt = simplePos.resConfs.filter { it.type == SimpleConfSpace.ResidueConf.Type.WildType }
			val compiledWt = compiledPos.confs.filter { it.id.startsWith("wt-") }
			if (simpleWt.size == compiledWt.size) {
				for ((rc, conf) in simpleWt.sortedBy { it.index }.zip(compiledWt.sortedBy { it.index })) {
					rcToConf[rc.index] = conf.index
				}
			} else {
				issues.add("wild-type count mismatch legacy=${simpleWt.size} compiled=${compiledWt.size}")
			}

			PositionMapping(simplePos, compiledPos, rcToConf, issues)
		}
	}

	private fun printMappings(mappings: List<PositionMapping>) {
		val complete = mappings.count { it.complete }
		println("mapping: completePositions=$complete/${mappings.size}")
		for (mapping in mappings) {
			val compiledName = mapping.compiledPos?.name ?: "<missing>"
			println("mapping.pos legacy=${mapping.simplePos.resNum} compiled=$compiledName complete=${mapping.complete} mapped=${mapping.rcToConf.size}/${mapping.simplePos.resConfs.size}")
			for (issue in mapping.issues) {
				println("  issue: $issue")
			}
			if (mapping.complete) {
				println("  rcToConf=${mapping.rcToConf.entries.joinToString(";") { "${it.key}->${it.value}" }}")
			}
		}
	}

		private fun printSampleAssignment(simple: SimpleConfSpace, mappings: List<PositionMapping>) {

			val input = System.getProperty("osprey.audit.input")?.takeIf { it.isNotBlank() } ?: return
		val designId = System.getProperty("osprey.bench.designId", "unknown")
		val state = System.getProperty("osprey.audit.state", "Complex")
		val assignments = loadFirstAssignments(input, designId, state) ?: run {
			println("sample.assignment: no matching row in $input")
			return
		}

		println("sample.assignment legacy=${assignments.joinToString(";")}")
		if (assignments.size != simple.positions.size) {
			println("sample.assignment issue: length=${assignments.size} expected=${simple.positions.size}")
			return
		}
		if (mappings.any { !it.complete }) {
			println("sample.assignment issue: cannot translate because mapping is incomplete")
			return
		}

		val compiledAssignments = IntArray(mappings.size)
		for (mapping in mappings) {
			val simpleRc = assignments[mapping.simplePos.index]
			val compiledConf = mapping.rcToConf[simpleRc] ?: run {
				println("sample.assignment issue: no mapping for ${mapping.simplePos.resNum} rc=$simpleRc")
				return
			}
			compiledAssignments[mapping.compiledPos!!.index] = compiledConf
		}
			println("sample.assignment compiled=${compiledAssignments.joinToString(";")}")
		}

		private fun runCudaSweep(
			simple: SimpleConfSpace,
			compiled: CompiledConfSpace,
			mappings: List<PositionMapping>
		) {
			if (mappings.any { !it.complete }) {
				println("cuda.sweep: skipped because mapping is incomplete")
				return
			}
			if (!CudaConfEnergyCalculator.isSupported()) {
				println("cuda.sweep: skipped because compiled CUDA is not supported on this node")
				return
			}

			val input = System.getProperty("osprey.audit.input")?.takeIf { it.isNotBlank() } ?: run {
				println("cuda.sweep: skipped because osprey.audit.input is not set")
				return
			}
			val designId = System.getProperty("osprey.bench.designId", "unknown")
			val state = System.getProperty("osprey.audit.state", "Complex")
			val maxConfs = Math.max(1, Integer.getInteger("osprey.compiledAudit.maxConfs", 512))
			val warmupConfs = Math.max(0, Integer.getInteger("osprey.compiledAudit.warmupConfs", 64))
			val repeats = Math.max(1, Integer.getInteger("osprey.compiledAudit.repeats", 2))
			val streamsList = parseIntList(System.getProperty("osprey.compiledAudit.streamsList"), listOf(8, 16, 32, 64, 128))
			val batchSizes = parseIntList(System.getProperty("osprey.compiledAudit.batchSizes"), listOf(256, 512, 1024, 2048, 4096))
			val precision = parsePrecision(System.getProperty("osprey.compiledAudit.precision", "Float32"))

			val legacyAssignments = loadAssignments(input, designId, state, maxConfs)
			val compiledAssignments = legacyAssignments.mapNotNull { translateAssignment(simple, mappings, it) }
			if (compiledAssignments.isEmpty()) {
				println("cuda.sweep: skipped because no assignments could be translated from $input")
				return
			}

			val jobs = compiledAssignments.map { conf ->
				ConfEnergyCalculator.MinimizationJob(conf, PosInterDist.all(compiled, conf))
			}
			val gpus = CudaConfEnergyCalculator.getGpusInfos()
			if (gpus.isEmpty()) {
				println("cuda.sweep: skipped because no CUDA GPUs were reported")
				return
			}

			println("cuda.sweep: precision=$precision assignments=${jobs.size}/${legacyAssignments.size} maxConfs=$maxConfs warmupConfs=$warmupConfs repeats=$repeats")
			println("cuda.sweep: streams=${streamsList.joinToString(",")} batchSizes=${batchSizes.joinToString(",")}")
			for (gpu in gpus) {
				println("cuda.gpu $gpu")
			}

			val gpu = gpus.first()
			for (streams in streamsList) {
				for (batchSize in batchSizes) {
					val gpuStreams = listOf(CudaConfEnergyCalculator.GpuStreams(gpu, streams))
					val setupStartNs = System.nanoTime()
					val ecalc = CudaConfEnergyCalculator(compiled, precision, gpuStreams, batchSize.toLong())
					val setupMs = nanosToMs(System.nanoTime() - setupStartNs)
					try {
						val warmupJobs = jobs.take(Math.min(warmupConfs, jobs.size))
						val warmupMs = if (warmupJobs.isNotEmpty()) {
							timeMinimize(ecalc, warmupJobs, batchSize)
						} else {
							0.0
						}

						val times = ArrayList<Double>()
						var checksum = Double.NaN
						for (repeat in 1..repeats) {
							val ms = timeMinimize(ecalc, jobs, batchSize)
							times.add(ms)
							checksum = energyChecksum(jobs)
							println(
								"cuda.repeat precision=$precision streams=$streams batchSize=$batchSize repeat=$repeat n=${jobs.size} ms=${"%.3f".format(Locale.US, ms)} checksum=${"%.6f".format(Locale.US, checksum)}"
							)
						}

						val bestMs = times.minOrNull() ?: Double.NaN
						val avgMs = times.average()
						val throughput = jobs.size / (bestMs / 1000.0)
						println(
							"cuda.result precision=$precision streams=$streams batchSize=$batchSize n=${jobs.size} repeats=$repeats setup_ms=${"%.3f".format(Locale.US, setupMs)} warmup_ms=${"%.3f".format(Locale.US, warmupMs)} best_ms=${"%.3f".format(Locale.US, bestMs)} avg_ms=${"%.3f".format(Locale.US, avgMs)} per_conf_ms=${"%.5f".format(Locale.US, bestMs / jobs.size)} throughput=${"%.2f".format(Locale.US, throughput)} checksum=${"%.6f".format(Locale.US, checksum)}"
						)
					} finally {
						ecalc.close()
					}
				}
			}
		}

		private fun timeMinimize(
			ecalc: CudaConfEnergyCalculator,
			jobs: List<ConfEnergyCalculator.MinimizationJob>,
			batchSize: Int
		): Double {
			val startNs = System.nanoTime()
			var offset = 0
			while (offset < jobs.size) {
				val end = Math.min(offset + batchSize, jobs.size)
				ecalc.minimizeEnergies(jobs.subList(offset, end))
				offset = end
			}
			return nanosToMs(System.nanoTime() - startNs)
		}

		private fun energyChecksum(jobs: List<ConfEnergyCalculator.MinimizationJob>) =
			jobs.fold(0.0) { sum, job -> sum + job.energy }

		private fun nanosToMs(ns: Long) = ns / 1_000_000.0

		private fun parseIntList(text: String?, defaults: List<Int>) =
			text
				?.split(",", ";", " ", "\t", "\n")
				?.mapNotNull { it.trim().takeIf { value -> value.isNotEmpty() }?.toInt() }
				?.filter { it > 0 }
				?.takeIf { it.isNotEmpty() }
				?: defaults

		private fun parsePrecision(text: String): Structs.Precision {
			return when (text.toLowerCase(Locale.US)) {
				"f32", "float", "float32" -> Structs.Precision.Float32
				"f64", "double", "float64" -> Structs.Precision.Float64
				else -> Structs.Precision.valueOf(text)
			}
		}

		private fun loadAssignments(input: String, designId: String, state: String, maxConfs: Int): List<IntArray> {
			val out = ArrayList<IntArray>()
			for (file in collectInputFiles(File(input))) {
				file.bufferedReader().useLines { lines ->
					val iterator = lines.iterator()
					if (!iterator.hasNext()) return@useLines
					val header = parseCsvLine(iterator.next())
					val columns = header.withIndex().associate { it.value to it.index }
					while (iterator.hasNext() && out.size < maxConfs) {
						val row = parseCsvLine(iterator.next())
						if (row.isEmpty()) continue
						fun col(name: String) = columns[name]?.let { row.getOrNull(it) } ?: ""
						if (designId != "unknown" && col("design_id") != designId) continue
						if (!sameState(col("state"), state)) continue
						out.add(parseAssignments(col("assignments")))
					}
				}
				if (out.size >= maxConfs) break
			}
			return out
		}

		private fun translateAssignment(
			simple: SimpleConfSpace,
			mappings: List<PositionMapping>,
			legacyAssignments: IntArray
		): IntArray? {
			if (legacyAssignments.size != simple.positions.size) {
				return null
			}
			val compiledAssignments = IntArray(mappings.size)
			for (mapping in mappings) {
				val simpleRc = legacyAssignments[mapping.simplePos.index]
				val compiledConf = mapping.rcToConf[simpleRc] ?: return null
				compiledAssignments[mapping.compiledPos!!.index] = compiledConf
			}
			return compiledAssignments
		}

		private fun loadFirstAssignments(input: String, designId: String, state: String): IntArray? {
			val files = collectInputFiles(File(input))
		for (file in files) {
			file.bufferedReader().useLines { lines ->
				val iterator = lines.iterator()
				if (!iterator.hasNext()) return@useLines
				val header = parseCsvLine(iterator.next())
				val columns = header.withIndex().associate { it.value to it.index }
				while (iterator.hasNext()) {
					val row = parseCsvLine(iterator.next())
					if (row.isEmpty()) continue
					fun col(name: String) = columns[name]?.let { row.getOrNull(it) } ?: ""
					if (designId != "unknown" && col("design_id") != designId) continue
					if (!sameState(col("state"), state)) continue
					return parseAssignments(col("assignments"))
				}
			}
		}
		return null
	}

	private fun collectInputFiles(input: File): List<File> {
		return when {
			input.isFile -> listOf(input)
			input.isDirectory -> input.listFiles()
				?.filter { it.isFile && it.name.endsWith(".csv") && !it.name.endsWith(".summary.csv") }
				?.sortedBy { it.name }
				?: emptyList()
			else -> emptyList()
		}
	}

	private fun parseAssignments(text: String) =
		text.split(";")
			.filter { it.isNotBlank() }
			.map { it.trim().toInt() }
			.toIntArray()

	private fun sameState(a: String, b: String) =
		a.toLowerCase(Locale.US) == b.toLowerCase(Locale.US)

	private fun describePosition(pos: DesignPosition): String {
		val sourceAtoms = pos.sourceAtoms
			.sortedBy { it.name }
			.joinToString(",") { it.name }
		val anchorGroups = pos.anchorGroups
			.joinToString(";") { group ->
				group.joinToString("+") { describeAnchor(it) }
			}
		return "posType=${pos.type} sourceAtoms=[$sourceAtoms] anchorGroups=[$anchorGroups]"
	}

	private fun describeAnchor(anchor: Anchor) =
		when (anchor) {
			is Anchor.Single -> "single(${anchor.a.name},${anchor.b.name},${anchor.c.name})"
			is Anchor.Double -> "double(${anchor.a.name},${anchor.b.name},${anchor.c.name},${anchor.d.name})"
		}

	private fun parseCsvLine(line: String): List<String> {
		val out = ArrayList<String>()
		val buf = StringBuilder()
		var quoted = false
		var i = 0
		while (i < line.length) {
			val ch = line[i]
			when {
				ch == '"' && quoted && i + 1 < line.length && line[i + 1] == '"' -> {
					buf.append('"')
					i++
				}
				ch == '"' -> quoted = !quoted
				ch == ',' && !quoted -> {
					out.add(buf.toString())
					buf.setLength(0)
				}
				else -> buf.append(ch)
			}
			i++
		}
		out.add(buf.toString())
		return out
	}
}
