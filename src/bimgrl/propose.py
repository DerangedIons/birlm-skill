"""Propose the blocks from the structure, the invariants, the stiffness and the compartments,
each block with the finding that produced it; then compare with gotranx's own ``suggest-blocks``
and list the judgment calls

The principle, whatever the model: a state takes a large step safely only if something covers
it, the base's exponential update, an implicit block or a sub-step block. Every other state is
advanced by forward Euler, which limits the step to 2|Re λ|/|λ|² of the stiffest mode it owns. The
structures named below (a Ca pool with its buffers, a Markov chain, V with its fast gates) are the
common instances in cardiac models, and the reason to group states so that an invariant stays
exact; they are not the list of what may be blocked. The generic rule covers what they leave, and
``propose.md`` shows, mode by mode, what covers each stiff mode and the step limit that remains.

Two rule sets, one per base scheme of ``block_rush_larsen``.

RL base (``--base rl``, the default; BIRL-M of the blocked-rush-larsen survey). The base
updates only the Hodgkin-Huxley gates exponentially and every other state by forward Euler, so
every non-gate state with a fast mode needs a block:
0. Generic rule, applied after the structural ones: every stiff mode owned by states that are
   still on forward Euler (V excepted, which the V block covers) gets an implicit block of those
   owners, with the whole chain or binding cluster of any member; owners that share a mode share
   a block, and an owner that shares a mode with an existing implicit block joins it.
1. Markov chains: as rule 1 of the GRL base.
2. One calcium block: every Ca-carrying state (the confirmed Ca pools, the total-Ca networks
   through them, their binding clusters) in one implicit block, merged with every implicit
   block that shares a state with it (the RyR chain). Not filtered by stiffness: a Ca pool
   outside a block is forward Euler and becomes the limit.
3. V block: as rule 4 of the GRL base.
4. No dyad block: it is a judgment call, ready to paste, to try in a tuning round.
5. Nesting: as rule 6 of the GRL base.
The file carries ``base_scheme = "rush_larsen"``, the gates as ``stiff_states`` and
``freeze_boundary_fluxes = true``.

GRL base (``--base grl``; BIMGRL, the paper's "Choice of blocks"):
1. Every Markov chain of an ion channel or of the release channel (a closed or exactly
   conserving chain that feeds V, or that the compartments name as the release channel) is an
   implicit block. A chain that does not feed V is listed, not blocked (under the RL base the
   generic rule blocks it when it owns a stiff mode).
2. Every binding cluster (a Ca pool with its buffer states) that owns a stiff mode is an
   implicit block; the non-stiff ones (cytosolic TnC, CaM) stay explicit.
3. The SR and the junctional pool, coupled by release, share one block, so that the release
   conserves Ca; if the junction has no buffer states it joins the SR block as a bare state.
4. V with the INa (and INaL) gates: a proportional sub-step block, rate trigger on V with
   atol 0.5 mV, h_min 0.001 ms.
5. The dyad: junctional Ca with its buffers, the release channel states, and the SR with its
   buffers: a binary sub-step block, h_min 0.005 ms, rate triggers on the junctional Ca
   (rtol 0.02) and on each open state of the release channel (atol 1e-3).
6. An implicit block that shares states with a sub-step block is placed entirely inside it.
"""

from __future__ import annotations

from pathlib import Path

from . import stiffness
from .workspace import Workspace, md_table

V_ATOL = 0.5
V_H_MIN = 0.001
DYAD_H_MIN = 0.005
DYAD_RTOL = 0.02
OPEN_ATOL = 1e-3
BASES = ("rl", "grl")


def _toml_list(xs):
    return "[" + ", ".join(f'"{x}"' for x in xs) + "]"


def propose(ws: Workspace, dt: float = 0.1, stiff_threshold: float = -20.0, nest_v_markov: bool = False, base: str = "rl") -> dict:
    if base not in BASES:
        raise SystemExit(f"--base must be one of {', '.join(BASES)}, got {base!r}")
    meta = ws.meta
    comps = meta.compartments
    structure = ws.read_json("structure.json")
    inv = ws.read_json("invariants.json")
    # the default file holds every mode at or below the default threshold, so a stricter threshold
    # filters it; a looser one needs the probe written by `bimgrl stiffness --threshold` at that value
    stiff_file = f"{stiffness.artifact(max(stiff_threshold, stiffness.THRESHOLD))}.json"
    stiff = ws.read_json(stiff_file) if (ws.out / stiff_file).exists() else None
    feeds_v = set(structure["feeds_v"]) | {structure["voltage"]}
    kind_of = {s["name"]: s["kind"] for s in structure["states"]}

    stiff_owners: dict[str, float] = {}
    if stiff:
        for m in stiff["modes"]:
            if m["min_lam_re"] <= stiff_threshold:
                for n in m["owners"]:
                    stiff_owners[n] = min(stiff_owners.get(n, 0.0), m["min_lam_re"])

    if base == "rl":
        from .stepper import ws_gates

        gates = list(meta.rl_gates or ws_gates(ws.model))
        explicit = {s["name"] for s in structure["states"]} - set(gates)
    else:
        gates = []
        explicit = {s["name"] for s in structure["states"] if s["update_grl"] == "FE"}
    modes = [m for m in stiff["modes"] if m["min_lam_re"] <= stiff_threshold] if stiff else []

    implicit: list[dict] = []
    other_chains: list[list[str]] = []
    notes: list[str] = []
    judgment: list[str] = []
    release = set(comps.release_channel)

    # 1. Markov chains
    for g in inv["invariants"]:
        if g["level"] != "subgroup" or g["structure"] != "markov":
            continue
        members = g["members"]
        if g["class"] == "improper":
            judgment.append(f"chain {', '.join(members)} is flagged as not conserving (see invariants.md); fix the model before blocking it")
            continue
        is_channel = bool(set(members) & feeds_v) or bool(set(members) & release)
        if not is_channel:
            other_chains.append(members)
            continue
        implicit.append(dict(states=members, why=f"Markov chain ({g['class']}); rule 1"))

    cluster_of: dict[str, dict] = {}
    sr, junc = comps.sr, comps.junction
    if base == "rl":
        implicit = _calcium_block(implicit, inv, comps, cluster_of, judgment)
    else:
        # 2. binding clusters that own a stiff mode
        clusters = [g for g in inv["invariants"] if g["level"] == "subgroup" and g["structure"] == "binding"]
        for g in clusters:
            members = g["members"]
            owners = [n for n in members if n in stiff_owners]
            if g["class"] == "improper":
                judgment.append(f"binding cluster {', '.join(members)} is flagged as not conserving; fix the model before blocking it")
                continue
            if owners:
                lam = min(stiff_owners[n] for n in owners)
                block = dict(states=list(members), why=f"pool with its buffers, stiff mode Re λ {lam:.4g} /ms owned by {', '.join(owners)}; rule 2")
                implicit.append(block)
                for n in members:
                    cluster_of[n] = block
            else:
                notes.append(f"binding cluster {', '.join(members)} owns no stiff mode (threshold {stiff_threshold:g} /ms): left explicit")
        if stiff is None:
            probe = "" if stiff_threshold <= stiffness.THRESHOLD else f" --threshold {stiff_threshold:g}"
            judgment.append(f"no {stiff_file}: every binding cluster was treated as non-stiff; run `bimgrl stiffness{probe}` first")

        # 3. SR + junction share a block
        if sr and junc:
            b_sr, b_j = cluster_of.get(sr), cluster_of.get(junc)
            if b_sr is not None and b_j is not None and b_sr is not b_j:
                merged = dict(states=b_sr["states"] + b_j["states"], why=f"{b_sr['why']} + {b_j['why']}; joined so that the release conserves Ca (rule 3)")
                implicit = [b for b in implicit if b is not b_sr and b is not b_j] + [merged]
                for n in merged["states"]:
                    cluster_of[n] = merged
            elif b_sr is not None and b_j is None:
                b_sr["states"] = b_sr["states"] + [junc]
                b_sr["why"] += f"; {junc} joined as a bare pool so that the release conserves Ca (rule 3)"
                cluster_of[junc] = b_sr
            elif b_j is not None and b_sr is None:
                b_j["states"] = b_j["states"] + [sr]
                b_j["why"] += f"; {sr} joined as a bare pool so that the release conserves Ca (rule 3)"
                cluster_of[sr] = b_j
        elif junc and not sr:
            judgment.append("a junctional pool but no SR was named in meta.toml [compartments]: check")

    substep: list[dict] = []
    # 4. V block
    v = meta.voltage
    na = [g for g in comps.na_gates if kind_of.get(g) == "gate"]
    if na:
        substep.append(dict(states=[v, *na], h_min=V_H_MIN, mode="proportional",
                            rate_triggers=[dict(state=v, atol=V_ATOL)], why="V with the INa/INaL gates, refined in proportion to |dV/dt| (rule 4)"))
    else:
        judgment.append("no INa gates named in meta.toml [compartments].na_gates: no V block proposed")

    # 5. dyad
    if junc:
        dyad = [junc]
        if junc in cluster_of:
            dyad = list(cluster_of[junc]["states"])
        chain_blocks = [b for b in implicit if set(b["states"]) & release]
        for b in chain_blocks:
            dyad += [s for s in b["states"] if s not in dyad]
        dyad += [s for s in release if s not in dyad]
        if sr and sr not in dyad:
            dyad.append(sr)
            if sr in cluster_of:
                dyad += [s for s in cluster_of[sr]["states"] if s not in dyad]
            judgment.append(f"the SR ({sr}) is in the dyad block by default; drop it if the Ca transient does not improve (the paper kept it where it helped)")
        triggers = [dict(state=junc, rtol=DYAD_RTOL)] + [dict(state=o, atol=OPEN_ATOL) for o in comps.release_open if o in dyad]
        if not comps.release_open:
            judgment.append("no open states of the release channel named: the dyad block triggers only on the junctional Ca")
        dyad_block = dict(states=dyad, h_min=DYAD_H_MIN, mode="binary", rate_triggers=triggers,
                          why="junctional Ca with its buffers, the release channel and the SR, refined while release is active (rule 5)")
        if base == "grl":
            substep.append(dyad_block)
        else:
            dyad_block["why"] = "the dyad: the whole calcium block and the release channel, refined while release is active"
            judgment.append(
                "no dyad sub-step block under the RL base: try without it first. In the survey it made T-World 5x slower "
                "than the suggested blocks and costs at every step in tissue, but it won by far in ten Tusscher 2004/2006, "
                "Grandi 2010 and ToR-ORd. To try it in a tuning round, add:\n\n```toml\n"
                + "\n".join(_substep_toml(dyad_block)) + "\n```")
    else:
        notes.append("no junctional/dyadic pool: no dyad block (as for models whose Ca enters the cytosol directly)")

    # 0. generic rule (RL base): whatever stiff mode the structural rules left on forward Euler
    if base == "rl":
        if stiff:
            _generic_blocks(implicit, substep, modes, explicit, inv, v, judgment)
        else:
            judgment.append(f"no {stiff_file}: the generic rule did not run, so a stiff state outside the chains and the calcium block "
                            "would stay on forward Euler unnoticed; run `bimgrl stiffness` first")
    blocked = set().union(*(b["states"] for b in implicit))
    for members in other_chains:
        if not set(members) & blocked:
            reason = "does not feed V" if base == "grl" else (
                "does not feed V and owns no stiff mode on forward Euler" if stiff else "does not feed V, and there is no stiffness file to judge it")
            notes.append(f"chain {', '.join(members)} {reason}: left explicit")

    # 6. nesting: an implicit block that shares states with a sub-step block goes inside it
    for b in implicit:
        for s in substep:
            shared = set(b["states"]) & set(s["states"])
            if shared and not set(b["states"]) <= set(s["states"]):
                s["states"] += [x for x in b["states"] if x not in s["states"]]
                s["why"] += f"; absorbed the implicit block {', '.join(b['states'])} (rule 6)"
    vblock = next((s for s in substep if s["mode"] == "proportional" and v in s["states"]), None)
    if nest_v_markov:
        if vblock:
            for b in implicit:
                if set(b["states"]) & feeds_v and not set(b["states"]) <= set(vblock["states"]):
                    vblock["states"] += b["states"]
                    vblock["why"] += f"; nested the voltage-dependent chain {', '.join(b['states'])}"
    nesting = _nesting_candidates(stiff, stiff_threshold, v, vblock, feeds_v, implicit)

    # what covers each stiff mode, and the forward Euler limit the blocks leave
    coverage, residual = _coverage(modes, implicit, substep, explicit)
    if base == "rl":
        for r in coverage:
            if r["uncovered_share"] < 0.5:
                continue
            line = (f"the stiff mode of {', '.join(r['owners'])} (Re λ {r['lam_re']:.4g} /ms, forward Euler limit {r['fe_limit_ms']:.3g} ms) is still "
                    f"mostly on forward Euler after these blocks, through {', '.join(r['uncovered'])}")
            if v in r["uncovered"]:
                line += f". {v} is covered only by a V sub-step block: name the fast gates that act on it in meta.toml [compartments].na_gates"
            else:
                line += ". No single state owns enough of it for the generic rule; put the listed states in one implicit block if the sweep stops at this limit"
            judgment.append(line)

    # gotranx's own suggestion for comparison
    suggested, suggest_msg = None, ""
    try:
        from gotranx.block_detection import DetectionSettings, blocks_to_toml, detect_blocks

        base_kw = dict(base_scheme="rush_larsen", freeze_boundary_fluxes=True) if base == "rl" else {}
        det = detect_blocks(ws.model.ode, DetectionSettings(dt=dt, membrane_potential=v, **base_kw))
        suggested = det
        (ws.out / f"blocks_suggested_dt{dt:g}.toml").write_text(blocks_to_toml(det))
    except Exception as e:  # noqa: BLE001
        suggest_msg = f"gotranx suggest-blocks failed: {type(e).__name__}: {e}"

    # every disagreement with gotranx is the user's call at Gate 6; the dyad is not compared because
    # gotranx does not detect it. Paths stay relative so the text is the same in every workspace
    if suggested:
        # compare which states are implicit, then how the states both make implicit are grouped
        mine_b = [set(b["states"]) for b in implicit]
        theirs_b = [set(b.states) for b in suggested.blocks.implicit]
        mine_s, theirs_s = set().union(*mine_b), set().union(*theirs_b)
        differ = []
        if theirs_s - mine_s:
            differ.append(f"it makes {', '.join(sorted(theirs_s - mine_s))} implicit and the rules leave them explicit")
        if mine_s - theirs_s:
            differ.append(f"it leaves {', '.join(sorted(mine_s - theirs_s))} explicit and the rules make them implicit")
        both = mine_s & theirs_s
        mine_g = {tuple(sorted(b & both)) for b in mine_b if b & both}
        theirs_g = {tuple(sorted(b & both)) for b in theirs_b if b & both}
        if mine_g != theirs_g:
            fmt = lambda gs: " ".join(f"[{', '.join(g)}]" for g in sorted(gs))  # noqa: E731
            differ.append(f"it groups the states both make implicit differently (rules {fmt(mine_g - theirs_g)}; gotranx {fmt(theirs_g - mine_g)})")
        mine_v = next((set(s["states"]) for s in substep if s["mode"] == "proportional" and v in s["states"]), set())
        theirs_v = next((set(b.states) for b in suggested.blocks.substep if v in b.states), set())
        if mine_v != theirs_v:
            only = "; ".join(f"only {who}: {', '.join(sorted(d))}" for who, d in (("gotranx", theirs_v - mine_v), ("the rules", mine_v - theirs_v)) if d)
            differ.append(f"the V block differs ({only})")
        judgment += [f"gotranx suggest-blocks disagrees: {d}" for d in differ]
        if differ:
            labels = "BIRL-M,BIRL-M-auto" if base == "rl" else "BIMGRL,BIMGRL-auto"
            judgment.append(f"to settle the gotranx disagreements on numbers, `bimgrl sweep <model_dir> --scheme {labels} "
                            f"--suggested <model_dir>/bimgrl/blocks_suggested_dt{dt:g}.toml` measures both block sets")

    stiff_states = None
    if base == "rl":
        in_implicit = set().union(*(set(b["states"]) for b in implicit))
        stiff_states = [g for g in gates if g not in in_implicit]
    toml = render_toml(implicit, substep, meta.id, round_=1, base=base, stiff_states=stiff_states)
    # validate
    try:
        from gotranx.blocks import parse_blocks, validate_blocks
        import tomllib

        validate_blocks(ws.model.ode, parse_blocks(tomllib.loads(toml), strict=True))
        valid = "valid (gotranx validate_blocks)"
    except Exception as e:  # noqa: BLE001
        valid = f"INVALID: {e}"
        judgment.append(valid)

    result = dict(base=base, stiff_states=stiff_states, implicit=implicit, substep=substep, notes=notes, judgment=judgment, valid=valid, toml=toml, dt=dt,
                  coverage=coverage, residual=residual, nesting=nesting, have_stiffness=stiff is not None,
                  suggested=_describe_detection(suggested) if suggested else None, suggest_msg=suggest_msg)
    return result


def _generic_blocks(implicit: list[dict], substep: list[dict], modes: list[dict], explicit: set, inv: dict, v: str, judgment: list[str]) -> None:
    """The generic rule of the RL base. A stiff mode whose owners are still on forward Euler limits
    the step whatever the states are called, so those owners get an implicit block. Owners that share
    a mode share a block, because the fast mode belongs to the group; the chain or binding cluster of
    a member comes whole, so that its sum stays exact; and an owner that shares a mode with an
    existing implicit block joins that block. V is left to the V sub-step block"""
    covered = set().union(*(b["states"] for b in implicit), *(s["states"] for s in substep))
    groups: list[dict] = []
    for m in modes:
        own = [n for n in m["owners"] if n in explicit and n not in covered and n != v]
        if not own:
            continue
        hit = [g for g in groups if set(g["states"]) & set(own)]
        g = dict(states=[], modes=[])
        for h in hit:
            g["states"] += [n for n in h["states"] if n not in g["states"]]
            g["modes"] += h["modes"]
        g["states"] += [n for n in own if n not in g["states"]]
        g["modes"].append(m)
        groups = [x for x in groups if not any(x is h for h in hit)] + [g]

    subgroups = [s for s in inv["invariants"] if s["level"] == "subgroup" and s["structure"] in ("markov", "binding")]
    for g in groups:
        whole, structured = [], set()
        for s in subgroups:
            if not set(s["members"]) & set(g["states"]):
                continue
            kind = "chain" if s["structure"] == "markov" else "binding cluster"
            if s["class"] == "improper":
                judgment.append(f"the {kind} {', '.join(s['members'])} has a stiff state but is flagged as not conserving (see invariants.md): "
                                "only the stiff state was blocked; fix the model before blocking the rest")
                continue
            g["states"] += [n for n in s["members"] if n not in g["states"] and n != v]
            structured |= set(s["members"])
            whole.append(f"the whole {kind} {', '.join(s['members'])} ({s['class']}), so that its sum stays exact")
        top = min(g["modes"], key=lambda m: m["min_lam_re"])
        desc = (f"stiff mode Re λ {top['min_lam_re']:.4g} /ms at t = {top['t_at_min']:g} ms (forward Euler limit {1e3 * top['fe_limit_ms']:.3g} µs), "
                f"owned by {', '.join(f'{n} {p:.2f}' for n, p in top['mean_participation'].items())}")
        if whole:
            desc += "; with " + "; ".join(whole)
        linked = set(g["states"]) | {n for m in g["modes"] for n in m["owners"]}
        targets = [b for b in implicit if set(b["states"]) & linked]
        in_implicit = set().union(*(b["states"] for b in implicit))
        new = [n for n in g["states"] if n not in in_implicit]
        if not new:
            continue
        if targets:
            targets[0]["states"] += new
            targets[0]["why"] += f"; joined by {', '.join(new)}, which shares a fast mode with this block ({desc}; generic rule)"
            if len(targets) > 1:
                judgment.append(f"{', '.join(new)} shares a fast mode with {len(targets)} implicit blocks and joined the first; "
                                "merge the blocks if the sweep shows an error at their boundary")
        else:
            caveat = "" if set(new) <= structured else (
                ". Not a pool with its buffers or a chain, so the implicit solve is not guaranteed stable: check against the reference")
            implicit.append(dict(states=new, generic=True, why=(
                f"{desc}; these states are not gates, so the RL base steps them by forward Euler (generic rule)" + caveat)))


def _coverage(modes: list[dict], implicit: list[dict], substep: list[dict], explicit: set) -> tuple[list[dict], dict | None]:
    """For each stiff mode, what advances each owner once the blocks are in place, and the share of
    the mode left on forward Euler. The share is an upper bound: it sums each state's largest
    participation over the checkpoints. The mode with at least half on forward Euler and the smallest
    limit is the step limit the blocks leave (the P1 of `stiffness`, after the blocks)"""
    where: dict[str, str] = {}
    for i, b in enumerate(implicit, 1):
        for n in b["states"]:
            where[n] = f"implicit block {i}"
    for s in substep:
        for n in s["states"]:
            where[n] = f"sub-step block ({s['mode']})"

    def how(n):
        return where.get(n) or ("forward Euler" if n in explicit else "exponential (base)")

    rows = []
    for m in modes:
        part = m.get("max_participation") or m["mean_participation"]
        left = {n: x for n, x in part.items() if how(n) == "forward Euler"}
        rows.append(dict(owners=m["owners"], lam_re=m["min_lam_re"], fe_limit_ms=m["fe_limit_ms"], t=m["t_at_min"],
                         cover={n: how(n) for n in m["owners"]}, uncovered_share=min(1.0, sum(left.values())),
                         uncovered=sorted(left, key=left.get, reverse=True)))
    limiting = [r for r in rows if r["uncovered_share"] >= 0.5]
    return rows, (min(limiting, key=lambda r: r["fe_limit_ms"]) if limiting else None)


def _nesting_candidates(stiff: dict | None, threshold: float, v: str, vblock: dict | None, feeds_v: set, implicit: list[dict]) -> list[dict]:
    """States outside the V block that the upstroke may under-resolve: they share a fast mode with V,
    or they act on V and are fast only at the checkpoints where V is. Evidence for the user, not a
    rule: nothing is nested unless a sweep shows the state limiting the step or leaving its range"""
    if not stiff or vblock is None:
        return []
    inside = set(vblock["states"])
    v_times: set[float] = set()
    shared: dict[str, float] = {}
    seen: dict[str, list[tuple[float, float]]] = {}
    for e in stiff["checkpoints"]:
        for m in e["modes"]:
            if m["lam_re"] > threshold:
                continue
            own = [n for n, _ in m["owners"]]
            if v in own:
                v_times.add(e["t"])
                for n in own:
                    if n not in inside:
                        shared[n] = min(shared.get(n, 0.0), m["lam_re"])
            for n in own:
                seen.setdefault(n, []).append((e["t"], m["lam_re"]))
    found = {n: f"shares a fast mode with {v} (Re λ {lam:.4g} /ms)" for n, lam in shared.items()}
    for n, hits in seen.items():
        ts = {t for t, _ in hits}
        if n in inside or n in found or n not in feeds_v or not v_times or not ts <= v_times:
            continue
        found[n] = f"acts on {v} and is fast only while {v} is (t = {min(ts):g} to {max(ts):g} ms, Re λ {min(lam for _, lam in hits):.4g} /ms)"
    out: list[dict] = []
    for n, why in found.items():
        block = next((b["states"] for b in implicit if n in b["states"]), [n])
        if not any(c["states"] == block for c in out):
            out.append(dict(states=block, why=why if block == [n] else f"{n} {why}; it sits in an implicit block, which nests whole"))
    return out


def _describe_detection(det) -> dict:
    return dict(
        implicit=[list(b.states) for b in det.blocks.implicit],
        substep=[dict(states=list(b.states), mode=b.mode, h_min=b.h_min,
                      rate_triggers=[dict(state=t.state, atol=t.atol, rtol=t.rtol) for t in b.rate_triggers]) for b in det.blocks.substep],
        reasons={",".join(k): v for k, v in det.reasons.items()},
    )


def _substep_toml(s: dict) -> list[str]:
    trig = ", ".join("{ " + ", ".join(f"{k} = {v!r}" if k != "state" else f'state = "{v}"' for k, v in t.items()) + " }" for t in s["rate_triggers"])
    return [f"# {s['why']}", "[[substep_blocks]]", f"states = {_toml_list(s['states'])}", f"h_min = {s['h_min']}", f'mode = "{s["mode"]}"',
            f"rate_triggers = [{trig}]"]


def render_toml(implicit, substep, model_id: str, round_: int, base: str = "grl", stiff_states=None) -> str:
    rules = "the RL-base rules (BIRL-M)" if base == "rl" else "the paper's GRL-base rules (BIMGRL)"
    lines = [f"# Blocks for {model_id} (round {round_}; see NOTES.md), proposed by bimgrl propose from", f"# {rules}; every block carries the finding that produced it.", ""]
    if base == "rl":
        lines += ["# The RL base: the Hodgkin-Huxley gates are updated exponentially, every other state by forward",
                  "# Euler unless it is in an implicit block", 'base_scheme = "rush_larsen"', f"stiff_states = {_toml_list(stiff_states or [])}",
                  "# Each implicit block takes the fluxes across its boundary explicitly", "freeze_boundary_fluxes = true", ""]
    for b in implicit:
        lines += [f"# {b['why']}", "[[implicit_blocks]]", f"states = {_toml_list(b['states'])}", ""]
    for s in substep:
        lines += _substep_toml(s) + [""]
    return "\n".join(lines)


def _calcium_block(implicit: list[dict], inv: dict, comps, cluster_of: dict, judgment: list[str]) -> list[dict]:
    """The RL-base calcium block: the Ca pools, the total-Ca networks and binding clusters through
    them, and every implicit block that shares a state with it, as one implicit block"""
    pools = [p for p in (comps.cytosol, comps.junction, comps.subsarcolemma, comps.sr) if p]
    if not pools:
        judgment.append("no Ca pool is named in meta.toml [compartments]: no calcium block, so under the RL base every Ca state runs on forward Euler")
        return implicit
    states, sources = list(pools), []

    def add(members):
        states.extend(m for m in members if m not in states)

    for g in inv["invariants"]:
        if not set(g["members"]) & set(pools):
            continue
        if g["level"] == "network":
            if g["class"] == "open":
                judgment.append(f"the network {g['name']} through the Ca pools does not conserve (class open, see invariants.md): its states were not added to the calcium block")
                continue
            add(g["members"])
            sources.append(f"{g['name']} ({g['class']})")
        elif g["level"] == "subgroup" and g["structure"] == "binding":
            if g["class"] == "improper":
                judgment.append(f"binding cluster {', '.join(g['members'])} is flagged as not conserving; fix the model before blocking it")
                continue
            add(g["members"])
            sources.append(f"binding cluster {', '.join(g['members'])}")
    merged = []
    while True:
        more = [b for b in implicit if b not in merged and set(b["states"]) & set(states)]
        if not more:
            break
        for b in more:
            add(b["states"])
        merged += more
    why = f"one calcium block: the Ca pools {', '.join(pools)}"
    if sources:
        why += f", with {'; '.join(sources)}"
    if merged:
        why += f"; absorbed {'; '.join(b['why'] for b in merged)}"
    why += " (RL base: a Ca state outside an implicit block runs on forward Euler and limits the step)"
    block = dict(states=states, why=why)
    for n in states:
        cluster_of[n] = block
    return [b for b in implicit if b not in merged] + [block]


def to_markdown(r: dict) -> str:
    base = r.get("base", "grl")
    out = ["# Proposed blocks", "", f"Base: {'Rush-Larsen (BIRL-M), frozen boundary fluxes' if base == 'rl' else 'generalized Rush-Larsen (BIMGRL)'}", "",
           f"Validation: {r['valid']}", ""]
    if base == "rl":
        out += [f"Stiff states (exponential under the base): {', '.join(r['stiff_states']) or 'none'}", ""]
    out += ["## Implicit blocks", ""] + ([f"- {', '.join(b['states'])}: {b['why']}" for b in r["implicit"]] or ["- none"])
    out += ["", "## Sub-step blocks", ""] + ([f"- {', '.join(s['states'])} ({s['mode']}, h_min {s['h_min']}, triggers {s['rate_triggers']}): {s['why']}" for s in r["substep"]] or ["- none"])
    out += ["", "## Stiff modes and what covers them", ""]
    if not r.get("have_stiffness"):
        out.append("No stiffness file: run `bimgrl stiffness` first. Without it nothing here says what the blocks leave on forward Euler.")
    elif not r["coverage"]:
        out.append("No stiff mode at this threshold.")
    else:
        out.append("A state takes a large step safely only if something covers it: the base's exponential update, an implicit block or a sub-step "
                   "block. A mode left mostly on forward Euler limits the step to the value in its row.")
        out += ["", md_table(
            ["owners", "min Re λ (/ms)", "FE limit (ms)", "t at min (ms)", "what advances each owner", "share left on forward Euler (upper bound)"],
            [[", ".join(c["owners"]), c["lam_re"], c["fe_limit_ms"], c["t"], "; ".join(f"{n}: {h}" for n, h in c["cover"].items()),
              f"{c['uncovered_share']:.2f}" + (" UNCOVERED" if c["uncovered_share"] >= 0.5 else "")] for c in r["coverage"]])]
        res = r["residual"]
        out += ["", "Predicted forward Euler step limit with these blocks: " + (
            f"{res['fe_limit_ms']:.4g} ms, from the mode of {', '.join(res['owners'])} (still on forward Euler: {', '.join(res['uncovered'])})."
            if res else "none; no stiff mode is left mostly on forward Euler.")]
    out += ["", "## Nesting candidates (not applied)", ""]
    out += [f"- {', '.join(c['states'])}: {c['why']}" for c in r.get("nesting", [])] or ["- none"]
    if r.get("nesting"):
        out += ["", "Nest one in the V block only when a sweep shows it limiting the step or leaving its range: every nested state is solved again "
                "on each sub-step. `--nest-v-markov` nests every implicit block that acts on V."]
    out += ["", "## Left explicit (by rule)", ""] + ([f"- {n}" for n in r["notes"]] or ["- nothing"])
    out += ["", "## Judgment calls for the user", ""] + ([f"- {j}".replace("\n", "\n  ") for j in r["judgment"]] or ["- none"])
    out += ["", f"## gotranx suggest-blocks at dt {r['dt']:g}, for comparison", ""]
    if r["suggested"]:
        s = r["suggested"]
        out += [f"- implicit: {b}" for b in s["implicit"]] or ["- implicit: none"]
        out += [f"- sub-step: {b['states']} ({b['mode']}, h_min {b['h_min']}, {b['rate_triggers']})" for b in s["substep"]] or ["- sub-step: none"]
        out += ["", "Reasons: " + "; ".join(f"{k}: {v}" for k, v in s["reasons"].items())]
        mine_i = {tuple(sorted(b["states"])) for b in r["implicit"]}
        theirs_i = {tuple(sorted(b)) for b in s["implicit"]}
        out += ["", "Differences: " + (
            "; ".join([f"only here: {list(b)}" for b in mine_i - theirs_i] + [f"only gotranx: {list(b)}" for b in theirs_i - mine_i]) or "none in the implicit blocks")]
        out.append("gotranx does not detect the dyad; its V block comes from the numeric fast-gate rule, which is the cross-check of the INa gate list.")
    else:
        out.append(r["suggest_msg"] or "not run")
    return "\n".join(out)


def run(ws: Workspace, dt: float = 0.1, stiff_threshold: float = -20.0, nest_v_markov: bool = False, write: bool = True, base: str = "rl") -> dict:
    r = propose(ws, dt, stiff_threshold, nest_v_markov, base)
    ws.write_json("propose.json", {k: v for k, v in r.items() if k != "toml"})
    ws.write_md("propose.md", to_markdown(r))
    target = ws.model_dir / "blocks.toml"
    path = ws.out / "blocks_proposed.toml"
    path.write_text(r["toml"])
    if write and not target.exists():
        target.write_text(r["toml"])
        print(f"wrote {target}")
    else:
        print(f"wrote {path} ({target} exists and was left alone)" if target.exists() else f"wrote {path}")
    return r
