##### File Organization functions for General Analyses


##################################################################
# Purpose: Creates a directory for a project in my standard format
# Arguments: $1 (name) -> name of the project to create and populate
# Arguments: $2 (dir) -> where you want the directory to be created
# Return: True or False
## Notes: Use the notes.txt file to record which scripts have been run on what. Only record successful runs and record the date.
##################################################################
generate_project_folder() 
{
    local name=$1
    local dir=$2
    current_dir=$(pwd)
    cd $dir
    mkdir $name
    cd $name
    mkdir scripts
    mkdir genome
    mkdir analyses
    mkdir results
    mkdir raw_data
    touch notes.txt
    echo "Project created, $(date)" >> notes.txt
    cd $current_dir
}

