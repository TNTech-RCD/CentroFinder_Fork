import os, subprocess

# Snakemake injects a global object named "snakemake"
fasta = snakemake.input["fasta"]
out_dat = snakemake.output["dat"]
logfile = snakemake.log[0]
trf_params = snakemake.params["trf_params"]

# Ensure directories exist
os.makedirs(os.path.dirname(output[0]), exist_ok=True)
os.makedirs(os.path.dirname(log[0]), exist_ok=True)

# Work directory for TRF is the results dir for this file
results_dir = os.path.dirname(output[0])

# Make input path relative to results_dir
input_rel = os.path.relpath(input[0], results_dir)

trf_command = f"trf {input_rel} {TRF_PARAM_STRING}"

# This code modified from https://stackoverflow.com/questions/45613881/what-would-be-an-elegant-way-of-preventing-snakemake-from-failing-upon-shell-r-e
try:
    # Run TRF; this will raise CalledProcessError on non-zero exit codes
    proc_output = subprocess.check_output(
        trf_command,
        shell=True,
        cwd=results_dir,
        stderr=subprocess.STDOUT,
    )

    # Log normal output
    with open(log[0], "wb") as lf:
        lf.write(proc_output)
        lf.write(b"\nTRF exit code: 0\n")

except subprocess.CalledProcessError as exc:
    # Log TRF output and exit code even on non-zero exit
    with open(log[0], "wb") as lf:
        if exc.output:
            lf.write(exc.output)
        lf.write(f"\nTRF exit code: {exc.returncode}\n".encode())

    # If TRF did not produce the expected .dat file, THEN treat as failure
    if not os.path.exists(output[0]):
        raise

# Final safety check: make sure the .dat file exists
if not os.path.exists(output[0]):
    raise Exception(
        f"TRF failed to produce expected output file: {output[0]}"
    )
