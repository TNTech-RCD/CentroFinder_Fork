import os, subprocess, shlex

# Snakemake injects a global object named "snakemake"
fasta = snakemake.input["fasta"]
out_dat = snakemake.output["dat"]
logfile = snakemake.log[0]
trf_param_string = snakemake.params.trf_param_string

# Ensure directories exist
os.makedirs(os.path.dirname(out_dat), exist_ok=True)
os.makedirs(os.path.dirname(logfile), exist_ok=True)

# Work directory for TRF is the results dir for this file
results_dir = os.path.dirname(out_dat)

# Make input path relative to results_dir
input_rel = os.path.relpath(fasta, results_dir)

# Build TRF command
trf_args = shlex.split(str(trf_param_string))
trf_cmd = ["trf", input_rel] + trf_args

trf_command = f"trf {shlex.quote(input_rel)} {trf_param_string}"

# This code modified from https://stackoverflow.com/questions/45613881/what-would-be-an-elegant-way-of-preventing-snakemake-from-failing-upon-shell-r-e
try:
    # Run TRF; this will raise CalledProcessError on non-zero exit codes
    proc_output = subprocess.check_output(
        trf_cmd,
        cwd=results_dir,
        stderr=subprocess.STDOUT,
    )

    # Log normal output
    with open(logfile, "wb") as lf:
        lf.write(proc_output)
        lf.write(b"\nTRF exit code: 0\n")

except subprocess.CalledProcessError as exc:
    # Log TRF output and exit code even on non-zero exit
    with open(logfile, "wb") as lf:
        if exc.output:
            lf.write(exc.output)
        lf.write(f"\nTRF exit code: {exc.returncode}\n".encode())

    # If TRF did not produce the expected .dat file, THEN treat as failure
    if not os.path.exists(out_dat):
        raise

# Final safety check: make sure the .dat file exists
if not os.path.exists(out_dat):
    raise Exception(
        f"TRF failed to produce expected output file: {out_dat}"
    )
