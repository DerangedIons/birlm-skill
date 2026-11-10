#!/bin/sh
# Preflight for the birlm skill: what is installed, and what language the project is in.
# POSIX sh on purpose: this has to run on a machine that has no Python yet.
# It only reports. The decisions are in references/environment.md.
#
#   sh scripts/env-check.sh [project_dir]      (default: the current directory)

project=${1:-.}
skill=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

first_line() { "$@" 2>/dev/null | head -n 1; }

row() { printf '| %s | %s | %s |\n' "$1" "$2" "$3"; }

tool() { # name, role, version command...
    name=$1 role=$2
    shift 2
    if path=$(command -v "$1" 2>/dev/null); then
        row "$name" "$path ($(first_line "$@"))" "$role"
    else
        row "$name" "missing" "$role"
    fi
}

echo "## Tools"
echo
echo "| tool | found | needed for |"
echo "|---|---|---|"

tool uv "host: builds the skill's Python environment (and can install Python itself)" uv --version

if py=$(command -v python3 2>/dev/null || command -v python 2>/dev/null); then
    ver=$("$py" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])' 2>/dev/null)
    ok=$("$py" -c 'import sys; print("ok" if sys.version_info >= (3, 11) else "too old, 3.11 needed")' 2>/dev/null)
    row python "$py ($ver, $ok)" "host, only if uv is missing"
else
    row python "missing" "host, only if uv is missing"
fi

# the order of compiler() in src/bimgrl/stepper.py
cc_found=missing
for c in "$CC" gcc cc clang; do
    [ -n "$c" ] || continue
    if p=$(command -v "$c" 2>/dev/null) && "$c" --version >/dev/null 2>&1; then
        cc_found="$p ($(first_line "$c" --version))"
        break
    fi
done
row "C compiler" "$cc_found" "sweep and time (without one: the generated-stepper path)"

tool julia "Julia on-ramp (scripts/jl), Julia target" julia --version

matlab=$(command -v matlab 2>/dev/null)
if [ -z "$matlab" ]; then
    # macOS and Linux installs are usually not on PATH
    for m in /Applications/MATLAB_R*.app/bin/matlab /usr/local/MATLAB/R*/bin/matlab; do
        [ -x "$m" ] && matlab=$m
    done
fi
row matlab "${matlab:-missing}" "MATLAB source or target (run with -batch)"

tool octave "stand-in for MATLAB when it cannot be run" octave --version

if [ -x "$skill/.venv/bin/gotranx" ]; then
    row "skill environment" "$skill/.venv ($(first_line "$skill/.venv/bin/gotranx" --version))" "every bimgrl command"
elif [ -d "$skill/.venv" ]; then
    row "skill environment" "$skill/.venv without gotranx (run bimgrl setup)" "every bimgrl command"
else
    row "skill environment" "not built (Setup in SKILL.md)" "every bimgrl command"
fi

echo
echo "## Project markers in $project"
echo
echo "| language | markers at the root | source files (two directories deep) |"
echo "|---|---|---|"

markers() { # files that exist at the project root, comma separated
    out=
    for f in "$@"; do
        [ -e "$project/$f" ] && out="${out:+$out, }$f"
    done
    printf '%s' "${out:--}"
}

count() { # number of files with any of the given extensions, two levels deep
    n=0
    for ext in "$@"; do
        k=$(cd "$project" 2>/dev/null && find . -maxdepth 3 -name "*.$ext" -not -path '*/.*' -not -path '*/node_modules/*' 2>/dev/null | wc -l)
        n=$((n + ${k:-0}))
    done
    printf '%s' "$n"
}

row Julia "$(markers Project.toml JuliaProject.toml Manifest.toml)" "$(count jl)"
row Python "$(markers pyproject.toml setup.py setup.cfg requirements.txt environment.yml)" "$(count py)"
row "C / C++ / CUDA" "$(markers CMakeLists.txt Makefile meson.build)" "$(count c h cpp hpp cc cu cuh)"
row MATLAB "$(markers startup.m Contents.m)$(ls "$project"/*.prj >/dev/null 2>&1 && printf ', *.prj')" "$(count m mlx)"
row "CellML / .ode" "-" "$(count cellml ode)"
