# Second-tier coverage: 4z80

The user chose `4z80_flex_p0` to accompany the already-running p1. Do not add
higher 4z80 tiers or select this system for the final long-run exploration.
All 38 systems now have at least two completed, running, or planned tiers;
this does not mean that two results per system have already completed.

Management job `12827710` submitted MARK* job `12827711_0` on compsci with
64 CPU threads, 192 GiB memory, account grisman, and a 72-hour limit. The
runner checks for AMD EPYC 9554, matching the current p1 anchor. This is a
coverage run; the time limit is not a completion forecast.

The serial lane is now p1 (`12812439_7`) → 1b6c_p2 (`12812439_41`) →
4z80_p0 (`12827711_0`) → 4wyq_p3 (`12827032_2`) → 1b6c_p4 (`12827323_0`).
Only the never-started 4wyq dependency changed. Running work was preserved,
and the existing eight-compsci/four-fennario concurrency caps remain in place.

No additional PACK* job is needed. GPU101 already produced all 39 p0 rows in
1412.23 seconds on four RTX PRO 6000 GPUs; 33 rows have all three components
Estimated. The fresh MARK* input preflight matches the 39 PACK* sequences,
residue lists, and structure checksum. Incomplete estimates must remain
visible in any subsequent numerical comparison.

`registry.json` records the actual submission and pairing. `designs.tsv`
contains the exact definition. The frozen execution package and audit live at
`/usr/xtmp/lz280/markstar_two_tier_20261006/prep_12827710`.
`restore_4z80.py` is a one-time management script and refuses to run again
while the submission registry exists.
