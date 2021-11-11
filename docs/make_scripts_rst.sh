# Loops through all the executable scripts and prints all their help text to scripts.rst so it can be used the the sphinx documentation


# Create the header of the main index page
echo "Welcome to vcstools's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :glob:

   modules
   scripts

The following scripts are provided as part of the vcstools package:
" > index.rst

echo ".. _scripts:

================
vcstools scripts
================

The following scripts are provided as part of the vcstools package:
" > scripts.rst

for vcstools_script in $(ls ../scripts/*py); do
    base_name=${vcstools_script##*/}
    echo $base_name
    header_dashes=$(perl -E "say '-' x ${#base_name}")

    # Split help into usage, decription and arguments
    python $vcstools_script -h | awk -v RS= '{print > ("help-" NR ".txt")}'
    usage=$(cat help-1.txt | awk '{print "  " $0}')
    description=$(cat help-2.txt | tr '\n' ' ')
    arguments=$(cat help-3.txt | awk '{print "  " $0}')
    rm help-*txt

    # Print the description to the main scripts page
    echo " - :ref:\`${base_name//_/-}-label\` - ${description}" >> index.rst

    # Print the help and pipe it to a file
    echo ".. _${base_name//_/-}-label:

${base_name}
${header_dashes}

${description}::


${usage}

${arguments}

"  >> scripts.rst
done