#!/usr/bin/bash

# Activate mamba mutatex environement
source /home/aliciapageau/mambaforge/etc/profile.d/conda.sh
source /home/aliciapageau/mambaforge/etc/profile.d/mamba.sh
mamba activate mutatex

# Deal with command line argument
VALID_ARGS=$(getopt -o iop:m: --long input,output,binding-interface,pdb-dir:,molecules-dir:,np: -- "$@")
if [[ $? != 0 ]]; then
    exit 1;
fi

eval set -- "$VALID_ARGS"

NP=1
INTERFACE=false
INPUT_DIR=$HOME/Documents/antifungal_project/mutateX
OUTPUT_DIR=`pwd` #usally run from output/mutatex

while true; do
  case "$1" in
    -i | --input ) INPUT_DIR="$2"; shift 2;;
    -o | --output ) OUTPUT_DIR="$2"; shift 2;;
    -p | --pdb-dir ) pdb="$2"; shift 2;;
    --binding-interface ) INTERFACE=true; shift;;
    -m | --molecules-dir ) jsonfiles="$2"; shift 2;;
    --np ) NP="$2"; shift 2;;
    -- ) shift; break;;
    * ) break ;;
  esac
done

for i in `ls $pdb | egrep '\.pdb$'` ; do
    echo "run mutateX scan for" $(basename $i)
    echo $NP
    
    # Prep directory tree with json molecules files
    mkdir -p $OUTPUT_DIR/mutations/$(basename $i .pdb)_model0_checked_Repair \
             $OUTPUT_DIR/repair/repair_$(basename $i .pdb)_model0_checked/molecules
    cp -a $jsonfiles/. $OUTPUT_DIR/repair/repair_$(basename $i .pdb)_model0_checked/molecules
    
    /usr/bin/bash $INPUT_DIR/seq3to1.sh $pdb/$i > $OUTPUT_DIR/$(basename $i | cut -d "_" -f1)_position_list.txt
    cd $OUTPUT_DIR/mutations/$(basename $i .pdb)_model0_checked_Repair
    xargs mkdir -p < $OUTPUT_DIR/$(basename $i .pdb)_position_list.txt
    for dir in *; do 
    	mkdir -- "$dir/molecules"
    	cp -a $jsonfiles/. $dir/molecules/
    done
    cd $OUTPUT_DIR
    

    if $INTERFACE; then
        echo 'Binding interface will be compute'
        mutatex $pdb/$i \
	--np $NP \
	-f suite5 \
	-R $INPUT_DIR/repair_runfile_template.txt \
	-M $INPUT_DIR/mutate_runfile_template.txt \
	-I $INPUT_DIR/interface_runfile_template.txt \
	--binding-energy \
	--foldx-log
    else
        echo 'Binding interface will not be compute'
        mutatex $pdb/$i \
	--np $NP \
	-f suite5 \
	-R $INPUT_DIR/repair_runfile_template.txt \
	-M $INPUT_DIR/mutate_runfile_template.txt \
	--foldx-log
    fi
    
    
    ddg2excel -p $(basename $i .pdb)_model0_checked.pdb \
    -d results/mutation_ddgs/$(basename $i .pdb)_model0_checked_Repair/ \
    -l $INPUT_DIR/mutation_list.txt \
    -F csv -o $(basename $i .pdb)_mutations_ddgs
    
    for binding in A-B A-C ; do
        ddg2excel -p $(basename $i .pdb)_model0_checked.pdb \
	-d results/interface_ddgs/$(basename $i .pdb)_model0_checked_Repair/${binding}/ \
	-l $INPUT_DIR/mutation_list.txt \
	-F csv -o $(basename $i .pdb)_${binding}_interfaces_ddgs
    done
done

mamba deactivate

cd $OUTPUT_DIR
python3 $INPUT_DIR/ddg_heatmap.py `ls *.csv`


#""""""""""""""""""""""""""""""""""""""""""""""""""""""""#
#Usage exemple :
#bash /home/aliciapageau/Documents/antifungal_project/00_script/new_mutatex.sh -p /home/aliciapageau/Documents/antifungal_project/top_candidates/ERG11/pdb/albicans/ -m /home/aliciapageau/Documents/antifungal_project/top_candidates/ERG11/output/mutatex/molecules --binding-interface
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""#
