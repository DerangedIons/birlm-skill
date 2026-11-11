# Asking the user: questions anyone on the team can answer

The people who run this skill know one side well and the other less: a cell modeller who has never chosen a time-integration scheme, a numerics or tissue-level person who has never needed to know what a dyad is. A gate question written in both vocabularies at once stops both of them. The scripts' reports stay precise and technical, because they are the record; what the agent *says* at a gate is a translation of them.

## The shape of every gate question

Write the question so that it can be answered by someone who knows neither field's jargon:

1. **What is being decided**, in one plain sentence. No state names, no scheme names, no symbols.
2. **What it changes**: accuracy, speed, or whether the result can be trusted. Give the size when a script has measured it.
3. **The default and why it is the default.** Put it first among the options, so that accepting it costs one keystroke.
4. **What happens if they do not know.** Always offer "not sure", and say what it does: take the default, or let a measurement decide. Wherever a sweep can compare the two choices (`BIRL-M` against `BIRL-M-auto`, a block with and without a member), offer that; the numbers then decide, not the user's background.
5. **Details last**, on a separate line for the reader who wants them: the state names, the eigenvalue, the rule number, the file.

Three rules follow from it:

- **No gate may require expertise the user did not claim.** If the only honest question needs a specialist ("is this state the open state of the release channel?"), ask the question behind it in terms of what can be seen ("this variable rises from 0 towards 1 while calcium is being released, and nothing else does: I will treat it as 'the release channel is open'. OK?") and keep "not sure, take your reading" as an option.
- **Say what a term means the first time it is used**, in half a sentence, from the glossaries below. Never assume that "stiff", "implicit", "dyad" or "upstroke" is known.
- **The same holds for what is reported**, not only for what is asked: the summary of the stiffness, the reason for each block, the hand-back paragraph. Say first what it means for the cell and for the step size, then name the modes and the states.

## Biology, for the reader who knows numerics

| term | what it is | what it means for the integration |
|---|---|---|
| membrane potential (V) | the voltage across the cell membrane | the variable every current depends on; it couples everything |
| upstroke | the first 1–2 ms of a beat, when V jumps by ~100 mV | the fastest event in the model, and brief: worth refining only while it lasts |
| fast sodium current (INa) and its gates | the current that drives the upstroke; its gates are on/off variables between 0 and 1 | the fastest gates; they have to be resolved together with V |
| gate (Hodgkin–Huxley) | a 0-to-1 variable relaxing towards a voltage-dependent target | a scalar linear ODE given V, so it has an exact exponential update: cheap and always stable |
| Markov chain of a channel | several states holding the probabilities of one channel's configurations; they sum to 1 | a linear system with a conserved sum; stepped one state at a time, the sum drifts |
| buffer | a molecule that binds calcium quickly and releases it again | a fast exchange between free and bound calcium: stiff, and their total must be kept |
| SR (sarcoplasmic reticulum) | the cell's internal calcium store | a large reservoir that empties into a tiny volume at every beat |
| junction, dyad, subspace | the tiny space the store releases into | small volume, so a very fast concentration; usually the stiffest calcium state |
| release channel (RyR) | the channel through which the store releases | switches on abruptly: a second brief fast event after the upstroke |
| crossbridge, troponin | the contraction machinery driven by calcium | often only an output of the model (nothing depends back on it), so an error there stays there |
| metabolites (ATP, NADH, the TCA-cycle species such as OAA and MAL) | the cell's energy chemistry | reaction networks with their own fast exchanges, unrelated to calcium or voltage |
| biomarker | one number read off a beat: peak voltage, beat duration (APD90), size of the calcium transient | the error measures; a scheme passes when they stay within a tolerance of the reference |

## Numerics, for the reader who knows the biology

| term | what it is | what it means for the cell model |
|---|---|---|
| time step (dt) | how far the simulation advances per update | larger is faster; too large and the result is wrong or blows up |
| stiff | some process in the model is much faster than what is being simulated | that process, not the physiology, dictates how small the step must be |
| forward Euler | the simplest update: new = old + step × rate | cheap, but unstable when the step is longer than the fastest process it handles |
| step limit | the largest step at which a method stays stable | the number the whole exercise is about raising |
| Rush–Larsen (RL) | an update that is exact for a gate | gates never limit the step; everything that is not a gate still can |
| mode, eigenvalue | a pattern of states that relaxes together, and how fast | tells which part of the cell is the fast process |
| owner of a mode | the states that carry most of that pattern | the states that need special treatment |
| implicit block | a group of states solved together, as a small linear system | removes the group's step limit and keeps its total (calcium, probability) exact; costs a small solve per step |
| sub-step block | a group of states that takes several small steps inside one large one, only while it is changing fast | resolves a brief fast event (the upstroke, the release) without slowing the rest of the beat |
| trigger tolerance | how much change per step is allowed before the sub-steps switch on | smaller means more accurate and slower |
| nesting | putting an implicit block inside a sub-step block | the group is solved again on every small step: safer during the fast event, slower |
| drift | a quantity that should stay constant (total calcium, probabilities summing to 1) slowly changing | the run looks fine and is wrong; the sign that a group was split |
| reference solution | the same model solved very accurately and slowly | what "correct" means when errors are measured |

## The gates, reworded

Examples of the shape, to adapt to the model at hand. The technical line of each comes from the script's report.

**Gate P (output language).** "Which language should the final, fast version of the model be written in? It changes nothing about the analysis, only the file you get at the end. Your project looks like Julia, so the default is Julia. *Details: gotranx generates C, C++, CUDA, Python, Julia and CytoZoo.*"

**Gate 0 (what is what in the model).** "I read the equations and guessed which variable plays which role in the cell; the table says what each role is. I am fairly sure of the rows without a mark. For the marked ones I show what I saw and what I concluded: tell me if any is wrong, or say 'not sure' and I keep my reading. A wrong guess here does not break the result: it makes the first set of blocks worse, and the test runs show it. *Details: the evidence column of `intake.md`.*" Then one short question per marked row, each of the form "this variable is the smallest calcium compartment and exchanges directly with the store, so I take it to be the space the store releases into. OK?"

**Gate 2 (something does not add up).** "Three variables of one channel should always sum to 1, and in the equations as written they do not: one transition is counted on one side only. This is in the model, not in the method. I can (a) point you at the two equations so you can check them against the paper [default], (b) go on as it is, in which case I leave this group alone and it may limit the step, or (c) stop. *Details: the flagged terms from `invariants.md`.*"

**Gate 6 (the blocks).** Present each block as what it is and what it buys, before its members: "The calcium in the cell's compartments and the molecules that bind it exchange much faster than a beat. Solving them together removes that limit and keeps total calcium exact. *Details: the states, the fastest rate, the step limit it removes.*" Then each judgment call as its own question:

- The dyad: "The calcium release is a second brief fast event. I can refine it with extra small steps, which is more accurate for the calcium transient and can be much slower, especially in tissue. Default: try without first, and add it only if the calcium error is too large."
- A disagreement with gotranx's own suggestion: "A second, automatic method chose a different grouping for five variables. Neither is known to be better for this model. Default: I run both and we keep the faster one that passes; it takes about N minutes."
- Nesting: "One channel reacts to voltage faster than the time step during the upstroke. Refining it together with voltage is safer but slower. Default: leave it as is, and revisit only if the test runs show it misbehaving."
- A generic block: "Two variables of the energy metabolism exchange very fast with each other. Left alone they would force a step about N times smaller for the whole model. I propose to solve the two together. This kind of group is less well-tried than the calcium one, so I check it against the accurate solution. OK?"
- A mode still uncovered: "After these blocks one fast process is still handled by the simple update, which caps the step at about N ms. Default: go on and see whether the test runs actually stop there."

**Gate 8 (a tuning round).** "At large steps the peak voltage is 3 mV too high: the upstroke is not resolved finely enough. I propose to let the voltage refinement switch on earlier, which costs a little time at each beat. Default: try it; I show the error and the cost before and after, and we keep it only if it helps. *Details: V-block trigger atol 0.5 → 0.1 mV.*"
