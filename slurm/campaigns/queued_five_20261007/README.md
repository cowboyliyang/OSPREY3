# Five approved MARK* additions

Management job `12832726` completed successfully on October 7, 2026. All five
MARK* submissions use account `grisman`, 64 CPUs, 192 GiB memory, and a 14-day
time limit. Their predecessor dependencies preserve the local concurrency
limits of nine compsci tasks and five fennario tasks.

| Design | Job | Location | Predecessor |
| --- | --- | --- | --- |
| `1gwc_p2` | `12832727_0` | fennario | `4wem_p10` (`12812439_2`) |
| `2xxm_p2` | `12832728_1` | fennario | `5dc4_p4` (`12812439_3`) |
| `3k3q_m1` | `12832729_2` | compsci | `3eb6_p6` (`12812439_13`) |
| `4znc_p5` | `12832730_3` | compsci | `1b6c_p2` (`12812439_41`) |
| `5d68_p4` | `12832731_4` | compsci | `5d68_p3` (`12826674_0`) |

The previously pending `4znc_p6` (`12812439_26`) now follows `4znc_p5`, so the
coverage chain is `1b6c_p2 -> 4znc_p5 -> 4znc_p6`. All 14 existing running
tasks retained their start times, commands, CPU counts, nodes, and time limits.
The five new jobs were verified PENDING after release from the submission
hold. The three conditional tiers remain deferred; DCC queues were not changed.

All five input preflights passed: structure checksums, positions, and complete
sequence sets match their measured PACK* cases. The frozen command templates
reuse the existing MARK* classpaths and record the expected CPU model: AMD
EPYC 9554 for compsci and Intel Xeon Gold 5320 for fennario. Older results retain
their original hardware records; a location change does not make them same-CPU
comparisons.

Preparation copied five structure files totaling 4,471,531 bytes. Runtime
artifacts and complete scheduler audit logs are under
`/usr/xtmp/lz280/markstar_queued_five_20261007/prep_12832726`; new run outputs use
the same scratch campaign root. Small versioned receipts are in
[`registry.json`](registry.json) and [`local_submissions.tsv`](local_submissions.tsv).
The application scripts reject duplicate application and do not restart
existing work.
