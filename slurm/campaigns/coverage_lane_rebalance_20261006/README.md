# Coverage lane rebalance

Applied by management job `12828045`. Prioritize the two coverage designs
behind the already-running `2hnu_p5` while retaining eight compsci and four
fennario serial lanes. Every changed task was pending with zero runtime.
The management job held those tasks during the five dependency changes,
verified the new dependencies, and released every temporary hold.

The three affected compsci lanes are now:

- `2hnu_p5 → 1b6c_p2 → 4znc_p6 → 1a0r_p0 → 2hnu_p11`
- `4z80_p1 → 4z80_p0 → 4wyq_p3 → 1b6c_p4`
- `3eb6_p6 → 3eb6_p7 → 4znc_p14`

All design names use `_flex_`. The runs use 64 CPU threads on EPYC 9554.
The two coverage designs no longer wait for `4z80_p1` and `3eb6_p6`.
`registry.json` records the applied changes; full before/after scheduler
records are under `/usr/xtmp/lz280/markstar_coverage_lane_rebalance_20261006/12828045`.

The user's schedule is two completed current-round tiers per system before
October 16, followed by a final 14-day round starting October 16. This change
advances two coverage tasks; it does not establish completion forecasts for
the full campaign or schedule the final round.
