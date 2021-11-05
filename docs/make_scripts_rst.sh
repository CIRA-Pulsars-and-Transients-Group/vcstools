# Loops through all the executable scripts and prints all their help text to scripts.rst so it can be used the the sphinx documentation

nl=$'\n'

# Create the header of the main scripts page
echo "# vcstools scripts${nl}" > scripts.md
echo "The following scripts are provided as part of the vcstools package:${nl}" >> scripts.md

for vcstools_script in $(ls ../scripts/*py); do
    base_name=${vcstools_script##*/}
    echo $base_name

    # Split help into usage, decription and arguments
    python $vcstools_script -h | awk -v RS= '{print > ("help-" NR ".txt")}'
    usage=$(cat help-1.txt)
    description=$(cat help-2.txt | tr '\n' ' ')
    arguments=$(cat help-3.txt)
    rm help-*txt

    # Print the description to the main scripts page
    echo " - [${base_name}](include_scripts/${base_name%%py}md) - ${description}" >> scripts.md

    # Print the help and pipe it to a file
    echo "# ${base_name}${nl}" > include_scripts/${base_name%%py}md
    echo "${description}${nl}"  >> include_scripts/${base_name%%py}md
    echo "\`\`\`"  >> include_scripts/${base_name%%py}md
    echo "${usage}${nl}"  >> include_scripts/${base_name%%py}md
    echo "${arguments}"  >> include_scripts/${base_name%%py}md
    echo "\`\`\`"  >> include_scripts/${base_name%%py}md
done