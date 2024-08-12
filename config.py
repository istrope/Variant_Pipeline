import os
import yaml
import shutil
import argparse

def check_tool(tool_name, user_path=None):
    """Check if a tool is in the PATH or if a user-provided path exists."""
    if user_path:
        if os.path.isfile(user_path):
            return user_path
        else:
            raise FileNotFoundError(f"Provided path for {tool_name} does not exist: {user_path}")
    else:
        tool_path = shutil.which(tool_name)
        if tool_path:
            return tool_path
        else:
            raise FileNotFoundError(f"{tool_name} is not installed and no path was provided.")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Generate a config.yaml file for Snakemake.")
    parser.add_argument("--bwa", help="Path to BWA executable.")
    parser.add_argument("--samtools", help="Path to Samtools executable.")
    parser.add_argument("--bcftools", help="Path to Bcftools executable.")
    parser.add_argument("--gatk", help="Path to GATK executable.")
    parser.add_argument("--varscan", help="Path to VarScan JAR file.", required=True)
    parser.add_argument("--reference_fasta", help="Path to reference FASTA file.", required=True)
    parser.add_argument("--funcotator_db", help="Path to Funcotator data sources.", required=True)
    parser.add_argument("--fastq_dir",help="Path to directory containing fastq files",required=True)
    parser.add_argument("--output", help="Output YAML file name.", default="config.yaml")
    parser.add_argument('--threads',help='number of threads for multi-threaded processes',default=4)

    args = parser.parse_args()

    # Create a dictionary to hold the paths
    config = {}

    # Check for each tool and add it to the config dictionary
    config["bwa"] = check_tool("bwa", args.bwa)
    config["samtools"] = check_tool("samtools", args.samtools)
    config["bcftools"] = check_tool("bcftools", args.bcftools)
    config["gatk"] = check_tool("gatk", args.gatk)
    config["varscan"] = check_tool("java", args.varscan)  # Assuming VarScan runs with java

    # Add reference fasta and funcotator db
    config["reference_fasta"] = args.reference_fasta
    config["funcotator_db"] = args.funcotator_db
    config['threads'] = args.threads
    config['fastq_dir']=args.fastq_dir

    # Write the config dictionary to a YAML file
    with open(args.output, "w") as yaml_file:
        yaml.dump(config, yaml_file, default_flow_style=False)

    print(f"{args.output} generated successfully.")

if __name__ == "__main__":
    main()
