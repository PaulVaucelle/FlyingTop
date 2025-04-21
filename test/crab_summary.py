import subprocess
import re

def get_crab_job_status(crab_directory):
    # Command to get crab status
    command = f"crab status -d {crab_directory}"
    
    try:
        # Execute the command and capture the output
        output = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        output = output.decode("utf-8")

        # Initialize job count variables
        job_counts = {
            'SUBMITTED': 0,
            'IDLE': 0,
            'RUNNING': 0,
            'COMPLETED': 0,
            'FAILED': 0
        }
                # Print the decoded output for inspection
        # print("Decoded output from 'crab status':")
        # print(output)

        # Extract job counts using regex (adjust this if the output format changes)
        submitted_match = re.search(r'submitted\s+(\d+)', output, re.IGNORECASE)
        idle_match = re.search(r'idle\s+(\d+)', output, re.IGNORECASE)
        running_match = re.search(r'running\s+(\d+)', output, re.IGNORECASE)
        completed_match = re.search(r'finished\s+(\d+)', output, re.IGNORECASE)
        failed_match = re.search(r'failed\s+(\d+)', output, re.IGNORECASE)

        if submitted_match:
            job_counts['SUBMITTED'] = int(submitted_match.group(1))
        if idle_match:
            job_counts['IDLE'] = int(idle_match.group(1))
        if running_match:
            job_counts['RUNNING'] = int(running_match.group(1))
        if completed_match:
            job_counts['COMPLETED'] = int(completed_match.group(1))
        if failed_match:
            job_counts['FAILED'] = int(failed_match.group(1))

        return job_counts

    except subprocess.CalledProcessError as e:
        print(f"Error running crab status: {e.output.decode('utf-8')}")
        return None

# Example usage
work_directories = [
    "crab_20240925_084230", "crab_20240925_084411", "crab_20240925_084501", 
    "crab_20240925_084547", "crab_20240925_084624", "crab_20240925_084800", 
    "crab_20240925_084840", "crab_20240925_084917", "crab_20240925_084951", 
    "crab_20240925_085032", "crab_20240925_085106", "crab_20240925_085200", 
    "crab_20240925_085234", "crab_20240925_085309", "crab_20240925_085343", 
    "crab_20240925_085419", "crab_20240925_085452", "crab_20240925_085524", 
    "crab_20240925_085554", "crab_20240925_085630", "crab_20240925_085704"
]


for work_directory in work_directories:
    crab_dir = f"./{work_directory}"
    job_status = get_crab_job_status(crab_dir)

    if job_status:
        print(f"CRAB Job Status for {work_directory}:")
        print(f"Submitted: {job_status['SUBMITTED']}")
        print(f"Idle: {job_status['IDLE']}")
        print(f"Running: {job_status['RUNNING']}")
        print(f"Completed: {job_status['COMPLETED']}")
        print(f"Failed: {job_status['FAILED']}")




