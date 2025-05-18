THREADS=16
SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
WORK_DIR=$(dirname ${SCRIPT_PATH})


while getopts "i:t:o:r" option; do
    case "${option}" in
        i) INPUT_FA=${OPTARG};;
        t) THREADS=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        r) RESTART=true;;
        *) exit 1;;
    esac
done

mkdir -p ${OUTPUT_DIR}/tmp/prodigal/checkpoint
DIR_FA=${INPUT_FA}
DIR_TMP=${OUTPUT_DIR}/tmp

# Create symbolic links only if they don't exist
for file in ${DIR_FA}/*fa; do
    symlink="${OUTPUT_DIR}/tmp/prodigal/$(basename ${file})"
    if [ ! -L "${symlink}" ]; then
        ln -s ${file} ${symlink}
    fi
done

# Check if all files in the directory end with .fa
invalid_files=false
for file in ${DIR_FA}/*; do
    if [ -f "$file" ] && [[ ! ${file} == *.fa ]]; then
        echo "Error: File ${file} does not end with .fa. Please rename the file."
        invalid_files=true
    fi
done

if [ "${invalid_files}" = true ]; then
    exit 1
fi

# Handle restart option
if [ "${RESTART}" = true ]; then
    rm -f ${OUTPUT_DIR}/tmp/prodigal/checkpoint/*.done
fi

# Process each file
for i in ${OUTPUT_DIR}/tmp/prodigal/*fa; do
    base=$(basename ${i} .fa)
    done_file="${OUTPUT_DIR}/tmp/prodigal/checkpoint/${base}.done"

    if [ -e "${done_file}" ]; then
        echo "${base} already processed."
        continue
    fi

    echo "Processing ${base}..."
    STARTTIME=$(date +%s)
    
    # Run prodigal
    prodigal -i ${i} \
             -o ${i%.fa}.gff \
             -a ${i%.fa}.faa \
             -d ${i%.fa}.fna \
             -f gff -p meta

    ENDTIME=$(date +%s)
    echo "[TIMER] Processing ${base} took $(($ENDTIME - $STARTTIME)) sec."

    # Create the .done file
    touch "${done_file}"
done

# -i /usr/commondata/public/gaoyang/yuanlin/test_for_NewVirFinder/experiment
# -o /usr/commondata/public/gaoyang/yuanlin/test_for_NewVirFinder/experiment/LRVM_results