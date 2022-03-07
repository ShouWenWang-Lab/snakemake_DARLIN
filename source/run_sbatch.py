import os
import argparse
# parse cmd line arguments
parser = argparse.ArgumentParser(description="run sbatch jobs")

# to pass a complex command (separated by space), please avoid the bash expansion by: --command '{command}'

parser.add_argument(
    "--job_name",
    type=str,
    default="name",
)

parser.add_argument(
    "--cores",
    type=str,
    default="name",
)

parser.add_argument(
    "--time",
    type=str,
    default="3",
    help="unit: hour",
)

parser.add_argument(
    "--mem",
    type=str,
    default="1G",
    help="requested memory",
)

parser.add_argument(
    "--command",
    type=str,
    default=".",
    help="the bash command to run",
)

args=parser.parse_args()
print("Run with sbatch jobs")
os.system("mkdir -p log")
sbatch_command=f"sbatch -p short -c {args.cores} -t {args.time}:00:00 --mem={args.mem} --job-name {args.job_name} --output=log/{args.job_name}-%j.o  --error=log/{args.job_name}-%j.e --mail-type=TIME_LIMIT_90,FAIL,END --wrap=\"{args.command}\""
print(f"submit job:   {sbatch_command}")
os.system(sbatch_command)