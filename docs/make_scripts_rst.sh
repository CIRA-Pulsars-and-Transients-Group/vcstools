# Loops through all the executable scripts and prints all their help text to an scripts.rst so it can be used the the sphinx documentation

nl=$'\n'

# Create the header of the main scripts page
echo "# vcstools scripts${nl}" > scripts.md
echo "The following scripts are provided as part of the vcstools package:${nl}" >> scripts.md

for vcstools_script in $(ls ../scripts/*py); do
    base_name=${vcstools_script##*/}
    echo $base_name

    # Print the help and pipe it to a file
    echo "# ${base_name}${nl}" > include_scripts/${base_name%%py}md
    python $vcstools_script -h >> include_scripts/${base_name%%py}md

    # Print the summary to the main scripts page
    echo " - [${base_name}](include_scripts/${base_name%%py}md) - " >> scripts.md
done