package osprey

import org.gradle.api.Project
import org.gradle.api.tasks.Exec
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.getValue
import java.io.File


fun Project.makeCudaTasks() {

	val compileCuda_forcefield by tasks.creating(Exec::class) {
		nvcc(this, "forcefield")
	}

	val compileCuda_ccd by tasks.creating(Exec::class) {
		nvcc(this, "ccd")
	}

	val compileCuda_residueForcefield by tasks.creating(Exec::class) {
		nvcc(this, "residueForcefield")
	}

	val compileCuda_residueCcd by tasks.creating(Exec::class) {
		nvcc(this, "residueCcd", maxRegisters=64)
	}

	val compileCuda_dp by tasks.creating(Exec::class) {
		nvcc(this, "dp")
	}

	val compileCuda_sampling by tasks.creating(Exec::class) {
		nvcc(this, "sampling")
	}

	@Suppress("UNUSED_VARIABLE")
	val compileCuda by tasks.creating {
		description = "Compile cuda kernels"
		dependsOn(
			compileCuda_forcefield,
			compileCuda_ccd,
			compileCuda_residueForcefield,
			compileCuda_residueCcd,
			compileCuda_dp,
			compileCuda_sampling
		)
	}
}

fun Project.nvcc(exec: Exec, kernelName: String, maxRegisters: Int? = null, profile: Boolean = false) {

	val nvcc = listOfNotNull(
		System.getenv("NVCC"),
		System.getenv("CUDA_HOME")?.let { "$it/bin/nvcc" },
		"/usr/local/cuda/bin/nvcc",
		"/usr/local/cuda-12.8/bin/nvcc",
		"nvcc"
	).firstOrNull { it == "nvcc" || File(it).canExecute() } ?: "nvcc"

	val args = mutableListOf(nvcc)

	if (profile) {
		// if profiling, compile for one arch with profiling/debug info
		// NOTE: change this to your GPU's arch
		args.addAll(listOf("-cubin", "-gencode=arch=compute_86,code=sm_86", "-lineinfo", "--ptxas-options=-v"))
	} else {
		// otherwise, compile for V100/Titan V and Ampere A5000 plus PTX for driver JIT compatibility
		args.addAll(listOf("-fatbin",
			"-gencode=arch=compute_70,code=sm_70",
			"-gencode=arch=compute_86,code=sm_86",
			"-gencode=arch=compute_86,code=compute_86"
		))
	}

	if (maxRegisters != null) {
		args.addAll(listOf("-maxrregcount", "$maxRegisters"))
	}

	args.addAll(listOf("$kernelName.cu", "-o", "$kernelName.bin"))

	exec.workingDir = file("src/main/resources/gpuKernels/cuda")
	exec.commandLine(args)
}
