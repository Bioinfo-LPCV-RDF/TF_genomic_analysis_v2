# ------------------------------------------------------------------------------
# ---------- prep_annotation
# ------------------------------------------------------------------------------
prepare_gff=$PATH_TO_COMPIL/bin/prepare_gff.py # DEPRECATED
generate_BedFromGff=$PATH_TO_COMPIL/bin/generate_bed_gff.py # DEPRECATED
bedops=/home/312.3-StrucDev/312.3.1-Commun/bedops/bin/bedops # v2.4.38 # DEPRECATED
bedmap=/home/312.3-StrucDev/312.3.1-Commun/bedops/bin/bedmap # v2.4.38 # DEPRECATED

# compute_space
spacing_mk=$PATH_TO_COMPIL/bin/get_interdistances.py
tffm_all_scores=/home/312.6-Flo_Re/312.6.1-Commun/scripts/get_all_score_tffm.py
spacing_mk_v2=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis/get_interdistances.py
plot_spacing_v2=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis/plot_spacings.r

compute_space(){

# compute_space -p <PEAKS> -ns <NEG_SET_1> -nb <NB_of_NS> -m <MATRICES> -n <NAMES> -od <RESULT_DIR> -th <THRESHOLDS> -max <MAXY>
local pocc=false
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks files set to: ${2}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "negative files set to: ${2}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "matrix files set to: ${2}";shift 2;;
		-n)
			local names=("${!2}")
			echo "names associated set to: ${2}";shift 2;;
        -th)
			local thresholds=("${!2}")
			echo "thresholds set to: ${2}";shift 2;;
        -maxy)
			local maxy=$2
			echo "maximum enrichment to display set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
        -g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -nb)
			local number_of_NS=$2
			echo "number of negative files to use set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-pc)
			local pocc=true
			echo "Pocc mode activated, pfm will be used to compute PWM score and Pocc";shift 1;;
		-h)
			usage compute_space; exit;;
		--help)
			usage compute_space; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_space; exit;;
	esac
done

local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $negative_sets ]; then echo "ERROR: -ns argument needed"; Errors+=1; fi
if [ -z $matrices ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $thresholds ]; then echo "ERROR: -th argument needed"; Errors+=1; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi

if [ -z $maxy ]; then echo "-maxy argument not used, no limits"; local maxy=0; fi
if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; local maxSpace=50; fi
if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; local minSpace=0; fi
if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; local offset_left=0; fi
if [ -z $offset_right ]; then echo "-or argument not used, no offset"; local offset_right=0; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi
if [ -z $number_of_NS ]; then echo "-nb argument not used, 1 negative set used"; local number_of_NS=1; fi
# if [ -z $colors ]; then echo "-od argument not used, "; local colors=('#40A5C7' '#F9626E' '#F0875A' '#307C95' '#BB4A52' '#B46544'); fi TODO implement

if [ $Errors -gt 0 ]; then usage compute_space; exit 1; fi




mkdir -p $results
# spacing_mk=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/get_interdistances.py

for ((i=0;i<${#peaks[@]};i++))
do
	local negative_set=${negative_sets[i]}
	local name=${names[i]} # TODO replace \s in name by _ 
	local matrice=${matrices[i]}
	local peak=${peaks[i]}
	
	local th=`join_by " " "${thresholds[@]}"`
	neg_set="NA"
	if [[ $peak != *".fa"* ]]; then
        bedtools getfasta -fi $genome -bed $peak -fo $results/pos_set.fa
        local peak=$results/pos_set.fa
    fi
	for ((j=1;j<$number_of_NS;j++))
	do
        if [[ ${negative_set} != *".fa"* ]];then
            bedtools getfasta -fi $genome -bed ${negative_set/_1_neg.bed/_${j}_neg.bed} -fo $results/neg_set_${j}.fa
            if [[ ${neg_set} == "NA" ]]; then
                neg_set="$results/neg_set_${j}.fa"
            else 
                neg_set+=" $results/neg_set_${j}.fa"
            fi
        else
			if [[ ${neg_set} == "NA" ]]; then
				neg_set="  ${negative_set/_1_neg.fa/_${j}_neg.fa}"
			else
            	neg_set+="  ${negative_set/_1_neg.fa/_${j}_neg.fa}"
			fi
        fi
        
	done
	echo $matrice
	if [[ $matrice == *".pfm"* ]]; then
        python $spacing_mk -mat $matrice -o $results -n ${4%.*}_pwm.svg -minInter $minSpace -maxInter $maxSpace  -pos $results/sets/interdist_peaks_pos.fas -th $th -neg $neg_set -points True -no_absolute_panel -ol $offset_left -or $offset_right --write_inter -one_panel -maxy $maxy
	fi
	if [[ $matrice == *".xml"* ]]; then
        python $spacing_mk -tffm $matrice -o $results -n ${name}_tffm.svg -minInter $minSpace -maxInter $maxSpace  -pos $peak -th $th -neg $neg_set -points True -no_absolute_panel -ol $offset_left -or $offset_right --write_inter -one_panel -maxy $maxy
    fi
	inkscape -z -e $results/${name}_tffm.png -w 1800 -h 1080 $results/${name}_tffm.svg
done



}


compute_space_v2(){

# compute_space -p <PEAKS> -ns <NEG_SET_1> -nb <NB_of_NS> -m <MATRICES> -n <NAMES> -od <RESULT_DIR> -th <THRESHOLDS> -maxy <MAXY> -maxs -mins -ol -or -g -wi -nap -op -co
local write_inter=false
local no_absolute_panel=false
local one_panel=false
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks files set to: ${2}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "negative files set to: ${2}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "matrix files set to: ${2}";shift 2;;
		-n)
			local names=("${!2}")
			echo "names associated set to: ${2}";shift 2;;
		-th)
			local thresholds=("${!2}")
			echo "thresholds set to: ${2}";shift 2;;
		-maxy)
			local maxy=$2
			echo "maximum enrichment to display set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
        -g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -nb)
			local number_of_NS=$2
			echo "number of negative files to use set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-wi)
			local write_inter=true
			echo "All spacings will be reported in file";shift 1;;
		-nap)
			local no_absolute_panel=true
			echo "absolute enrichment panel won't be plotted";shift 1;;
		-op)
			local one_panel=true
			echo "relative enrichment panels will be merged";shift 1;;
		-co)
			local colors=("${!2}")
			echo "Colors to use: ${2}"; shift 2;;
		-h)
			usage compute_space_v2; exit;;
		--help)
			usage compute_space_v2; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_space; exit;;
	esac
done

local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $negative_sets ]; then echo "ERROR: -ns argument needed"; Errors+=1; fi
if [ -z $matrices ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
if [ -z $names ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $thresholds ]; then echo "ERROR: -th argument needed"; Errors+=1; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi

if [ -z $maxy ]; then echo "-maxy argument not used, no limits"; local maxy=0; fi
if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; local maxSpace=50; fi
if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; local minSpace=0; fi
if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; local offset_left=0; fi
if [ -z $offset_right ]; then echo "-or argument not used, no offset"; local offset_right=0; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi
if [ -z $number_of_NS ]; then echo "-nb argument not used, 1 negative set used"; local number_of_NS=1; fi
if [ -z $colors ]; then echo "-co argument not used, "; local colors=('#004899' '#004088' '#003876' '#0050aa' '#1961b2' '#3272bb' '#4c84c3' '#6696cc' '#7fa7d4' '#99b9dd'); fi

if [ $Errors -gt 0 ]; then usage compute_space_v2; exit 1; fi


local additional_python=""
local additional_R=()
if [ $write_inter = true ];then
	additional_python+=("-wi")
fi
if [ $no_absolute_panel = true ];then
	additional_python+=("-nap")
	additional_R+=("--no_absolute_panel")
fi
if [ $one_panel = true ];then
	additional_R+=("--one_panel")
fi
local list_additional_R=`join_by " " "${additional_R[@]}"`
local list_additional_python=`join_by " " "${additional_python[@]}"`
local list_colors=`join_by "," "${colors[@]}"`

for ((i=0;i<${#peaks[@]};i++))
do
	
	local negative_set=${negative_sets[i]}
	local name=${names[i]// /_}
	local matrice=${matrices[i]}
	local peak=${peaks[i]}
	local Neg_is_fasta="false"
	local length_mat=0
	mkdir -p $results/${name}/scores
	
	local th=`join_by " " "${thresholds[@]}"`
	neg_set="NA"
	if [[ $peak == *"score"* ]];then
		local pos_file=$peak # no preparation needed
	elif [[ $peak != *".fa"* ]]; then
		if [ ! -f $results/${name}/${name}_pos_set.fa ];then
			bedtools getfasta -fi $genome -bed $peak -fo $results/${name}/${name}_pos_set.fa # need to tranform in fasta
		fi
		local peak=$results/${name}/${name}_pos_set.fa
	fi
	if [[ $peak == *".fa"* ]]; then # need to compute scores
		if [[ $matrice == *".xml"* ]]; then
			if [ ! -f ${results}/${name}/scores/${name}_tffm_scores_pos.tsv ];then
				python $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_pos.tsv -pos $peak -t ${matrice}
			fi
			local pos_file=${results}/${name}/scores/${name}_tffm_scores_pos.tsv
		elif [[ $matrice == *".pfm"* ]]; then
			local length_mat=$(awk -v OFS="[\t ]" '{if(NR==1){if($5=="SIMPLE"){typM="simple"} else {typM="dependency"};count=-1;next};if(typM=="simple"){count++};if(typM=="dependency"){if($1=="DEPENDENCY"){count-=2;exit};count++}}END{print count}' $matrice)
			
			if [ ! -f ${results}/${name}/scores/$(basename $peak).scores ];then
				python $scores_prog -m ${matrice} -f $peak -o ${results}/${name}/scores/
			fi
			local pos_file=${results}/${name}/scores/$(basename $peak).scores
		fi
	fi

	for ((j=1;j<=$number_of_NS;j++))
	do
		local tmp_fasta="NA"
		if [[ ${negative_set} == *"score"* ]];then
			local tmp_neg_set="${negative_set/_1_neg.fa.scores/_${j}_neg.fa.scores}"
		elif [[ ${negative_set} == *".bed"* ]];then
			if [ ! -f $results/${name}/${name}_set_${j}_neg.fa ];then
				bedtools getfasta -fi $genome -bed ${negative_set/_1_neg.bed/_${j}_neg.bed} -fo $results/${name}/${name}_set_${j}_neg.fa
			fi
			local tmp_fasta=$results/${name}/${name}_set_${j}_neg.fa
			local Neg_is_fasta="true"
		elif [[ ${negative_set} == *".fa"* ]];then
			local Neg_is_fasta="true"
			local tmp_fasta=${negative_set/_1_neg.fa/_${j}_neg.fa}
		fi
		if [ $Neg_is_fasta == "true" ]; then
			if [[ $matrice == *".xml"* ]]; then
				if [ ! -f ${results}/${name}/scores/${name}_tffm_scores_${j}_neg.tsv ];then
					python $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_${j}_neg.tsv -pos $tmp_fasta -t ${matrice}
				fi
				local tmp_neg_set="${results}/${name}/scores/${name}_tffm_scores_${j}_neg.tsv"
			elif [[ $matrice == *".pfm"* ]]; then
				if [ ! -f ${results}/${name}/scores/$(basename $tmp_fasta).scores ];then
					python $scores_prog -m ${matrice} -f $tmp_fasta -o ${results}/${name}/scores/
				fi
				local tmp_neg_set="${results}/${name}/scores/$(basename $tmp_fasta).scores"
			fi
		fi
		if [[ ${neg_set} == "NA" ]]; then
			local neg_set="$tmp_neg_set"
		else
			local neg_set+=" $tmp_neg_set" 
		fi
	done 
 	python $spacing_mk_v2 -o $results/${name}/${name} -smax $maxSpace -smin $minSpace -pos $pos_file -th $th -neg $neg_set -ol $offset_left -or $offset_right -lm $length_mat $list_additional_python
#  	python $spacing_mk_v2 -o $results/${name} -smax $maxSpace -smin $minSpace -pos $pos_file -th $th -neg $neg_set -ol $offset_left -or $offset_right -lm $length_mat $list_additional_python ###### Modif Laura
	/home/prog/R/R-3-5-0/bin/Rscript $plot_spacing_v2 -e ${results}/${name}/${name}_enrichment.tsv -d ${results}/${name}/${name} -r ${results}/${name}/${name}_rate.tsv -m $minSpace -c $list_colors -y $maxy $list_additional_R
done

}


# v17/11/2021
compute_space(){

# compute_space -p <PEAKS> -ns <NEG_SET_1> -nb <NB_of_NS> -m <MATRICES> -n <NAMES> -od <RESULT_DIR> -th <THRESHOLDS> -max <MAXY>
local pocc=false
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks files set to: ${peaks[@]}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "negative files set to: ${negative_sets[@]}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "matrix files set to: ${matrices[@]}";shift 2;;
		-n)
			local names=("${!2}")
			echo "names associated set to: ${names[@]}";shift 2;;
        -th)
			local thresholds=("${!2}")
			echo "thresholds set to: ${thresholds[@]}";shift 2;;
        -maxy)
			local maxy=$2
			echo "maximum enrichment to display set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
        -g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -nb)
			local number_of_NS=$2
			echo "number of negative files to use set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-pc)
			local pocc=true
			echo "Pocc mode activated, pfm will be used to compute PWM score and Pocc";shift 1;;
		-h)
			usage compute_space; exit;;
		--help)
			usage compute_space; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_space; exit;;
	esac
done

local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $negative_sets ]; then echo "ERROR: -ns argument needed"; Errors+=1; fi
if [ -z $matrices ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
if [ -z $name ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $thresholds ]; then echo "ERROR: -th argument needed"; Errors+=1; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi

if [ -z $maxy ]; then echo "-maxy argument not used, no limits"; local maxy=0; fi
if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; local maxSpace=50; fi
if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; local minSpace=0; fi
if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; local offset_left=0; fi
if [ -z $offset_right ]; then echo "-or argument not used, no offset"; local offset_right=0; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi
if [ -z $number_of_NS ]; then echo "-nb argument not used, 1 negative set used"; local number_of_NS=1; fi
# if [ -z $colors ]; then echo "-od argument not used, "; local colors=('#40A5C7' '#F9626E' '#F0875A' '#307C95' '#BB4A52' '#B46544'); fi TODO implement

if [ $Errors -gt 0 ]; then usage compute_space; exit 1; fi




mkdir -p $results
# spacing_mk=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/get_interdistances.py

for ((i=0;i<${#peaks[@]};i++))
do
	local negative_set=${negative_sets[i]}
	local name=${names[i]} # TODO replace \s in name by _ 
	local matrice=${matrices[i]}
	local peak=${peaks[i]}
	
	local th=`join_by " " "${thresholds[@]}"`
	neg_set="NA"
	if [[ $peak != *".fa"* ]]; then
        bedtools getfasta -fi $genome -bed $peak -fo $results/pos_set.fa
        local peak=$results/pos_set.fa
    fi
	for ((j=1;j<$number_of_NS;j++))
	do
        if [[ ${negative_set} != *".fa"* ]];then
            bedtools getfasta -fi $genome -bed ${negative_set/_1_neg.bed/_${j}_neg.bed} -fo $results/neg_set_${j}.fa
            if [[ ${neg_set} == "NA" ]]; then
                neg_set="$results/neg_set_${j}.fa"
            else 
                neg_set+=" $results/neg_set_${j}.fa"
            fi
        else
			if [[ ${neg_set} == "NA" ]]; then
				neg_set="  ${negative_set/_1_neg.fa/_${j}_neg.fa}"
			else
            	neg_set+="  ${negative_set/_1_neg.fa/_${j}_neg.fa}"
			fi
        fi
        
	done
	echo $matrice
	if [[ $matrice == *".pfm"* ]]; then
        python $spacing_mk -mat $matrice -o $results -n ${4%.*}_pwm.svg -minInter $minSpace -maxInter $maxSpace  -pos $results/sets/interdist_peaks_pos.fas -th $th -neg $neg_set -points True -no_absolute_panel -ol $offset_left -or $offset_right --write_inter -one_panel -maxy $maxy
	fi
	if [[ $matrice == *".xml"* ]]; then
        python $spacing_mk -tffm $matrice -o $results -n ${name}_tffm.svg -minInter $minSpace -maxInter $maxSpace  -pos $peak -th $th -neg $neg_set -points True -no_absolute_panel -ol $offset_left -or $offset_right --write_inter -one_panel -maxy $maxy
    fi
	inkscape -z -e $results/${name}_tffm.png -w 1800 -h 1080 $results/${name}_tffm.svg
done



}

# v17/11/2021
compute_space_v2(){

# compute_space -p <PEAKS> -ns <NEG_SET_1> -nb <NB_of_NS> -m <MATRICES> -n <NAMES> -od <RESULT_DIR> -th <THRESHOLDS> -maxy <MAXY> -maxs -mins -ol -or -g -wi -nap -op -co
local write_inter=false
local no_absolute_panel=false
local one_panel=false
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-p)
			local peaks=("${!2}")
			echo "peaks files set to: ${peaks[@]}";shift 2;;
		-ns)
			local negative_sets=("${!2}")
			echo "negative files set to: ${negative_sets[@]}";shift 2;;
		-m)
			local matrices=("${!2}")
			echo "matrix files set to: ${matrices[@]}";shift 2;;
		-n)
			local names=("${!2}")
			echo "names associated set to: ${names[@]}";shift 2;;
		-th)
			local thresholds=("${!2}")
			echo "thresholds set to: ${thresholds[@]}";shift 2;;
		-maxy)
			local maxy=$2
			echo "maximum enrichment to display set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
        -g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -nb)
			local number_of_NS=$2
			echo "number of negative files to use set to: ${2}";shift 2;;
		-od)
			local results=$2
			echo "output directory set to: ${2}";shift 2;;
		-mat_type)
			local matrix_type=$2
			echo "should be 'ASYMMETRIC' or 'SYMMETRIC': ${2}"; shift 2;;
		-wi)
			local write_inter=true
			echo "All spacings will be reported in file";shift 1;;
		-nap)
			local no_absolute_panel=true
			echo "absolute enrichment panel won't be plotted";shift 1;;
		-op)
			local one_panel=true
			echo "relative enrichment panels will be merged";shift 1;;
		-co)
			local colors=("${!2}")
			echo "Colors to use: ${colors[@]}"; shift 2;;
		-h)
			usage compute_space_v2; exit;;
		--help)
			usage compute_space_v2; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage compute_space; exit;;
	esac
done

local Errors=0
if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
if [ -z $negative_sets ]; then echo "ERROR: -ns argument needed"; Errors+=1; fi
if [ -z $matrices ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
if [ -z $names ]; then echo "ERROR: -n argument needed"; Errors+=1; fi
if [ -z $thresholds ]; then echo "ERROR: -th argument needed"; Errors+=1; fi
if [ -z $results ]; then echo "ERROR: -od argument needed"; Errors+=1; fi

if [ -z $maxy ]; then echo "-maxy argument not used, no limits"; local maxy=0; fi
if [ -z $maxSpace ]; then echo "-maxs argument not used, using 50"; local maxSpace=50; fi
if [ -z $minSpace ]; then echo "-mins argument not used, using 0"; local minSpace=0; fi
if [ -z $offset_left ]; then echo "-ol argument not used, no offset"; local offset_left=0; fi
if [ -z $offset_right ]; then echo "-or argument not used, no offset"; local offset_right=0; fi
if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas" ; fi
if [ -z $number_of_NS ]; then echo "-nb argument not used, 1 negative set used"; local number_of_NS=1; fi
if [ -z $colors ]; then echo "-co argument not used, "; local colors=('#004899' '#004088' '#003876' '#0050aa' '#1961b2' '#3272bb' '#4c84c3' '#6696cc' '#7fa7d4' '#99b9dd'); fi

if [ $Errors -gt 0 ]; then usage compute_space_v2; exit 1; fi


local additional_python=""
local additional_R=()
if [ $write_inter = true ];then
	additional_python+=("-wi")
fi
if [ $no_absolute_panel = true ];then
	additional_python+=("-nap")
	additional_R+=("--no_absolute_panel")
fi
if [ $one_panel = true ];then
	additional_R+=("--one_panel")
fi
local list_additional_R=`join_by " " "${additional_R[@]}"`
local list_additional_python=`join_by " " "${additional_python[@]}"`
local list_colors=`join_by "," "${colors[@]}"`

for ((i=0;i<${#peaks[@]};i++))
do
	
	local negative_set=${negative_sets[i]}
	local name=${names[i]// /_}
	local matrice=${matrices[i]}
	local peak=${peaks[i]}
	local Neg_is_fasta="false"
	local length_mat=0
	mkdir -p $results/${name}/scores
	
	local th=`join_by " " "${thresholds[@]}"`
	neg_set="NA"
	if [[ $peak == *"score"* ]];then
		local pos_file=$peak # no preparation needed
	elif [[ $peak != *".fa"* ]]; then
		if [ ! -f $results/${name}/${name}_pos_set.fa ];then
			bedtools getfasta -fi $genome -bed $peak -fo $results/${name}/${name}_pos_set.fa # need to tranform in fasta
		fi
		local peak=$results/${name}/${name}_pos_set.fa
	fi
	if [[ $peak == *".fa"* ]]; then # need to compute scores
		if [[ $matrice == *".xml"* ]]; then
			if [ ! -f ${results}/${name}/scores/${name}_tffm_scores_pos.tsv ];then
				python $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_pos.tsv -pos $peak -t ${matrice}
			fi
			local pos_file=${results}/${name}/scores/${name}_tffm_scores_pos.tsv
		elif [[ $matrice == *".pfm"* ]]; then
			local length_mat=$(awk -v OFS="[\t ]" '{if(NR==1){if($5=="SIMPLE"){typM="simple"} else {typM="dependency"};count=-1;next};if(typM=="simple"){count++};if(typM=="dependency"){if($1=="DEPENDENCY"){count-=2;exit};count++}}END{print count}' $matrice)
			
			if [ ! -f ${results}/${name}/scores/$(basename $peak).scores ];then
				python $scores_prog -m ${matrice} -f $peak -o ${results}/${name}/scores/
			fi
			local pos_file=${results}/${name}/scores/$(basename $peak).scores
		fi
	fi

	for ((j=1;j<=$number_of_NS;j++))
	do
		local tmp_fasta="NA"
		if [[ ${negative_set} == *"score"* ]];then
			local tmp_neg_set="${negative_set/_1_neg.fa.scores/_${j}_neg.fa.scores}"
		elif [[ ${negative_set} == *".bed"* ]];then
			if [ ! -f $results/${name}/${name}_set_${j}_neg.fa ];then
				bedtools getfasta -fi $genome -bed ${negative_set/_1_neg.bed/_${j}_neg.bed} -fo $results/${name}/${name}_set_${j}_neg.fa
			fi
			local tmp_fasta=$results/${name}/${name}_set_${j}_neg.fa
			local Neg_is_fasta="true"
		elif [[ ${negative_set} == *".fa"* ]];then
			local Neg_is_fasta="true"
			local tmp_fasta=${negative_set/_1_neg.fa/_${j}_neg.fa}
		fi
		if [ $Neg_is_fasta == "true" ]; then
			if [[ $matrice == *".xml"* ]]; then
				if [ ! -f ${results}/${name}/scores/${name}_tffm_scores_${j}_neg.tsv ];then
					python $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_${j}_neg.tsv -pos $tmp_fasta -t ${matrice}
				fi
				local tmp_neg_set="${results}/${name}/scores/${name}_tffm_scores_${j}_neg.tsv"
			elif [[ $matrice == *".pfm"* ]]; then
				if [ ! -f ${results}/${name}/scores/$(basename $tmp_fasta).scores ];then
					python $scores_prog -m ${matrice} -f $tmp_fasta -o ${results}/${name}/scores/
				fi
				local tmp_neg_set="${results}/${name}/scores/$(basename $tmp_fasta).scores"
			fi
		fi
		if [[ ${neg_set} == "NA" ]]; then
			local neg_set="$tmp_neg_set"
		else
			local neg_set+=" $tmp_neg_set" 
		fi
	done 
 	python $spacing_mk_v2 -o $results/${name}/${name} -smax $maxSpace -smin $minSpace -pos $pos_file -th $th -neg $neg_set -ol $offset_left -or $offset_right -lm $length_mat $list_additional_python
	#exit 0
#  	python $spacing_mk_v2 -o $results/${name} -smax $maxSpace -smin $minSpace -pos $pos_file -th $th -neg $neg_set -ol $offset_left -or $offset_right -lm $length_mat $list_additional_python ###### Modif Laura
	/home/prog/R/R-3-5-0/bin/Rscript $plot_spacing_v2 -e ${results}/${name}/${name}_enrichment.tsv -d ${results}/${name}/${name} -r ${results}/${name}/${name}_rate.tsv -m $minSpace -c $list_colors -y $maxy $list_additional_R
	  sed "1d" ${results}/${name}/${name}_spacing_pos.tsv | sed "s/-/:/" | sed "s/:/\t/g" | sed "s/_/\t/" | awk -v OFS="\t" '{print $0,$3-$2+1}' | sed "1ichr\tStart\tEnd\tConf\tSpace\tScore1\tScore2\tmatricePosition1\tmatricePosition2\tcorrectedPosition1\tcorrectedPosition2\tSize" > ${results}/${name}/${name}_spacing_pos_reform.tsv
	  head ${results}/${name}/${name}_spacing_pos_reform.tsv
	#ajout Zscore spacing 25 aout 2021 + modif nov 2021
	echo $matrix_type
	echo $th
	#echo $(which R) # defaultis /home/312.3-StrucDev/312.3.1-Commun/bin/sl6/R -> /home/prog/R/R-3.2.5/bin/R but I need 4.1
	/home/prog/R/R-4.1.0/bin/Rscript $Zscore_spacing ${results}/${name}/${name}_spacing_pos_reform.tsv $matrix_type $th ${results}/${name}
done

}

# usage of deprecated functions
case
	compute_space) 
echo -e "
==========
usage: compute_space -p <LIST of FILE> -ns <LIST of FILE> -m <LIST of FILE> 
       -n <LIST of STRING> -od <RESULT_DIR> -th <LIST of FLOAT> 
       -nb <INT> -maxy [INT] -maxs [INT] -mins [INT] -ol [INT] -or [INT]
       -g [FILE] [-pc]

general infos: computes spacing analysis on peak files and negatives sets 
               using PWM/TFFM.

-- Mandatory arguments:
    -p     FILE     :    List of peak files (bed, narrowpeak format). Example:
                         (\"FILE1\" \"FILE2\" ... \"FILEn\")
    -ns    FILE     :    List of names of the negative sets. Example:
                         (\"FILE1\" \"FILE2\" ... \"FILEn\")
    -nb    INT      :    number of negatives sets to use against each peak 
                         file (suggested: 3; be sure to create enough negative
                         file using compute_NS)
    -m     FILE     :    List of matrices files. Example:
                         (\"FILE1\" \"FILE2\" ... \"FILEn\")
    -n     STRING   :    List of names for the sets in the ROCs. Example:
                         (\"FILE1\" \"FILE2\" ... \"FILEn\")
    -od    PATH     :    Set the Output directory
    -th    FLOAT    :    list of thresholds to use. (between 0-1 for TFFM
                         matrices & between -60-0 for PWMs). Choosing 
                         thresholds can be quite complicated. One way to get
                         satisfying threshold is to set them empirically and
                         adjust them according to the results. If dots are 
                         missing for some spacing, thresholds are too high, 
                         reduce them (closer to 0 for TFFM; closer to -60 for 
                         PWM). If enrichment increase and reduce itself 
                         rapidly, everywhere in the graph, you are seeing 
                         noises, increase you thresholds (closer to 1 for 
                         TFFM; closer to 0 for PWM). Example:
                         (\"FLOAT1\" \"FLOAT2\" ... \"FLOATn\")

-- optional arguments :
    -maxy  INT      :    maximum enrichment to display on y-axis.(default: NA)
    -maxs  INT      :    maximum spacing to compute. For dimers (using 
                         monomeric matrices), a value between 30 to 50 bp 
                         seems a good start. For tetramers (using dimeric 
                         matrices) a value between 70 to 100 bp is a good 
                         start. You may adjust those values after your first
                         results. (default: )
    -mins  INT      :    minimum spacing to compute. Usually set to 0. You may
                         change this value if you add offsets, this will
                         prevent huge drop off in enrichment for the first 
                         spacings computed.
    -ol    INT      :    offset to apply on the left of the matrix used. This
                         argument (in addition to -or) allows to change the 
                         way spacings are counted. This may be usefull when 
                         you want to count space around a consensus sequence
                         or from the center of the matrix.
    -or    INT      :    offset to apply on the right of the matrix used. See
                         -ol option for details.
    -g     FILE     :    Fasta of the genome used as reference.
    -pc             :    pocc mode.
\n";;
	compute_space_v2) # compute_space -p <PEAKS> -ns <NEG_SET_1> -nb <NB_of_NS> -m <MATRICES> -n <NAMES> -od <RESULT_DIR> -th <THRESHOLDS> -maxy <MAXY> -maxs -mins -ol -or -g -wi -nap -op -co
echo -e "
==========
usage: compute_space -p <LIST of FILE> -ns <LIST of FILE> -m <LIST of FILE> 
       -n <LIST of STRING> -od <RESULT_DIR> -th <LIST of INT> 
       -nb <NB_of_NS> -maxy [INT] -maxs [INT] -mins [INT] -ol [INT] -or [INT]
       -g [FILE] -co [LIST of STRING] [-nap] [-op] [-wi]

general infos: computes spacing version 2 analysis on peak files and negatives
               sets using PWM/TFFM. This new version used R to plot the 
               results. Plots are refreshed & codes are consolidated.

-- Mandatory arguments:
    -p     FILE     :    List of peak files (bed, narrowpeak format,fasta or     
                         scores files accepted). Example:
                         (\"FILE1\" \"FILE2\" ... \"FILEn\")
    -ns    FILE     :    List of names of the negative set (format to use for 
                         name: *_1_neg.*). Example:
                         (\"FILE1\" \"FILE2\" ... \"FILEn\")
    -nb    INT      :    number of negatives sets to use against each peak 
                         file (suggested: 3; be sure to create enough negative
                         file using compute_NS)
    -m     FILE     :    List of matrices  (either PWM or TFFM). Example:
                         (\"FILE1\" \"FILE2\" ... \"FILEn\")
    -n     STRING   :    List of names for the sets in the outputs. Example:
                         (\"STRING1\" \"STRING2\" ... \"STRINGn\")
    -od    PATH     :    Set the output directory
    -th    FLOAT    :    list of thresholds to use. (between 0-1 for TFFM
                         matrices & between -60-0 for PWMs). Choosing 
                         thresholds can be quite complicated. One way to get
                         satisfying thresholds is to set them empirically and
                         adjust them according to the results. A warning will 
                         be printed in terminal if thresholds are too high, 
                         reduce them (closer to 0 for TFFM; closer to -60 for 
                         PWM). If enrichment increase and reduce itself 
                         rapidly, everywhere in the graph, you are seeing 
                         noises, increase you thresholds (closer to 1 for 
                         TFFM; closer to 0 for PWM). Example:
                         (\"FLOAT1\" \"FLOAT2\" ... \"FLOATn\")

-- optional arguments :
    -maxy  INT      :    maximum enrichment to display on y-axis.(default: NA)
    -maxs  INT      :    maximum spacing to compute. For dimers (using 
                         monomeric matrices), a value between 30 to 50 bp 
                         seems a good start. For tetramers (using dimeric 
                         matrices) a value between 70 to 100 bp is a good 
                         start. You may adjust those values after your first
                         results. (default: 50)
    -mins  INT      :    minimum spacing to compute. Usually set to 0. You may
                         change this value if you add offsets, this will
                         prevent huge drop off in enrichment for the first 
                         spacings computed.
    -ol    INT      :    offset to apply on the left of the matrix used. This
                         argument (in addition to -or) allows to change the 
                         way spacings are counted. This may be useful when 
                         you want to count space around a consensus sequence
                         or from the center of the matrix.
    -or    INT      :    offset to apply on the right of the matrix used. See
                         -ol option for details.
    -g     FILE     :    Fasta of the genome used as reference.
    -op             :    activate \"one panel\" mode. ER, DR and IR for both 
                         negative & positive sets will be add to each other 
                         and plotted in only one panel. Thi mode can be useful
                         if your matrix is palindromic.
    -nap            :    Activate \"no absolute panel\" mode. Absolute 
                         Enrichment will not be plotted on the side of the 
                         panel(s). Can be called with \"-op\" argument.
    -wi             :    Activate \"write interdistance\" mode. This will 
                         trigger the report of every conformation & spacing 
                         found for each file treated. Reports only duos of 
                         sites with score higher than the lowest threshold 
                         used.
    -co    COLORs   :    Colors attributed to each threshold, in same order. 
                         Colors accepted are Hexadecimals (default blu).
\n";;

esac



cooking_meth_bk_21juin2021(){

while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-g)
			local genome=$2
			echo "Fasta of the genome set to: ${2}";shift 2;;
        -p)
			local peaksCov=$2
			echo "Tab-seaparated file for bound region with chr, start, end, covDAP, covAMP: ${2}";shift 2;;
        -s)
			local pfmResults=$2
			echo "pfm search results set to: ${2}";shift 2;;
		-l)
			local motifLength=$2
			echo "Motif length set to: ${2}";shift 2;;
		-c)
			local cutoff=$2
			echo "PFM score cutoff set to: ${2}";shift 2;;
		-c2)
			local cutoff2=$2
			echo "second step PFM score cutoff set to: ${2}";shift 2;;
		-sym)
			local symmetry=$2
			echo "[yes] or [no] to specifies wether or not the motif is symmetric: ${2}";shift 2;;
		-m)
			local matrix_pfm=$2
			echo "TF pfm matrix: ${2}";shift 2;;
        -o)
			local outdir=$2
			echo "output directory set to: ${2}";shift 2;;
		-h)
			usage cooking_meth; exit;;
		--help)
			usage cooking_meth; exit;;
		*)
			echo "Error in arguments"
			echo $1; usage cooking_meth; exit;;
	esac
done

local Errors=0
if [ -z $peaksCov ]; then echo "ERROR: -p argument needed";Errors+=1;fi
if [ -z $genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
if [ -z $pfmResults ]; then echo "ERROR: -s argument needed";Errors+=1;fi
if [ -z $motifLength ]; then echo "ERROR: -l argument needed";Errors+=1;fi
if [ -z $cutoff ]; then echo "ERROR: -c argument needed";Errors+=1;fi
if [ -z $cutoff2 ]; then echo "ERROR: -c2 argument needed";Errors+=1;fi
if [ -z $outdir ]; then echo "ERROR: -o argument needed";Errors+=1;fi
if [ -z $symmetry ]; then echo "ERROR: -sym argument needed";Errors+=1;fi
if [ $Errors -gt 0 ]; then echo "error somewhere"; exit 1; fi

time $full_methylation -p $peaksCov -g $genome -m $methMap -s $pfmResults -o $outdir -c $cutoff -l $motifLength

# the Rscript below will:
# compute methylation statistics
# plot these statistics
# and save these statistics as vectors into an RData object, so that it can be load into an R env to create pretty figures !

$plot_meth_full $outdir $motifLength $symmetry $cutoff2

$figs_meth_violin $outdir $matrix_pfm $symmetry 
}


heatmap_reads(){ # DEPRECATED
local centered="False"
# heatmap_reads
while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-b)
			local bedgraphs=("${!2}")
			echo "bedgraphs files set to: ${bedgraphs[@]}";shift 2;;
        -n)
			local names=("${!2}")
			echo "names set to: ${names[@]}";shift 2;;
		-od)
			local out_dir=$2
			echo "output directory set to: ${2}";shift 2;;
		-or)
			local order=$2
			echo "order set to: ${2}";shift 2;;
		-s)
			local size_file=$2
			echo "file of chromosome sizes set to: ${2}";shift 2;;
		-p)
			local peaks=$2
			echo "peaks set to: ${2}";shift 2;;
		-ws)
			local window_size=`calc $2/2`
			echo "window size around center of peaks set to: ${2}";shift 2;;
		-ow)
			local order_window=$2
			echo "window for ordering: $2"; shift 2;;
		-c)
			local centered="True";
			echo "plots will be centered"; shift 1;;
		-h)
			usage heatmap_reads; return;;
		--help)
			usage heatmap_reads; return;;
		*)
			echo "Error in arguments"
			echo $1; usage heatmap_reads; return;;
	esac
done

local Errors=0
if [ -z $bedgraphs ]; then echo "ERROR: -b argument needed or need to be a LIST";Errors+=1;fi
if [ -z $names ]; then echo "ERROR: -n argument needed or need to be a LIST";Errors+=1;fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi
if [ -z $size_file ]; then echo "ERROR: -s argument needed";Errors+=1;fi
if [ -z $peaks ]; then echo "ERROR: -p argument needed";Errors+=1;fi

if [ ${#names[@]} -ne ${#bedgraphs[@]} ]; then
	echo "ERROR: -n & -b have to be lists of same length"; Errors+=1
fi

if [ -z $window_size ]; then echo "-ws argument not used, using 1000bp";local window_size=500 ;fi
if [ -z $order ]; then echo "-or argument not used, ordering by CFR";local order=3 ;fi
if [ -z $order_window ]; then echo "-ow argument not used, ordering on full sequences";local order_window="0:$window_size" ;fi
if [[ $order_window != *":"* ]]; then
	echo ""
	Errors+=1
fi

if [ $Errors -gt 0 ]; then usage heatmap_reads; return 1; fi

mkdir -p $out_dir

if [ $window_size -gt 0 ]; then 
awk -v len=$window_size -v OFS="\t" '{print $1, int($2+($3-$2)/2-len), int($2+($3-$2)/2+len),$4}' ${peaks} > $out_dir/${names[0]}_${names[1]}_peaks.bed
else
cp ${peaks} $out_dir/${names[0]}_${names[1]}_peaks.bed
fi
for ((i=0;i<${#bedgraphs[@]};i++))
do
	if [ ! -f $out_dir/${names[i]}.bw ]; then
		$bdg_to_bw ${bedgraphs[i]} $size_file $out_dir/${names[i]}.bw
	fi
done
processed=$out_dir/${names[0]}_${names[1]}_peaks.bed
local names=`join_by " " "${names[@]}"`

$Python_TFFM $heatmap_mk -b $out_dir -p $processed -r $out_dir -s $order -n $names -w $window_size -c $centered --orderwindow $order_window


}


spacing_impact(){ # DEPRECATED

while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n)
			local names=("${!2}")
			echo "names set to: ${names[@]}";shift 2;;
		-cp)
			local comp_dir=${2}
			echo "initial comparison directory set to: ${2}";shift 2;;
		-od)
			local out_dir=${2}
			echo "output directory set to: ${2}";shift 2;;
		-g)
			local genome=${2}
			echo "Genome fasta set to: ${2}";shift 2;;
		-m)
			local matrix=${2}
			echo "matrix set to: ${2}";shift 2;;
        -maxs)
			local maxSpace=$2
			echo "maximum spacing to compute set to: ${2}";shift 2;;
        -mins)
			local minSpace=$2
			echo "minimum spacing to compute set to: ${2}";shift 2;;
        -ol)
			local offset_left=$2
			echo "offset on the left set to: ${2}";shift 2;;
        -or)
			local offset_right=$2
			echo "offset on the right set to: ${2}";shift 2;;
		-sth)
			local starting_th=$2
			echo "starting threshold set to: ${2}";shift 2;;
		-eth)
			local ending_th=$2
			echo "ending threshold set to: ${2}";shift 2;;
		-ith)
			local increment_th=$2
			echo "threshold increments set to: ${2}";shift 2;;
		-sp)
			local spacing=("${!2}")
			echo "spacings set to: ${spacing[@]}";shift 2;;
		-c)
			local colors=("${!2}")
			echo "colors set to: ${colors[@]}";shift 2;;
		-h)
			usage spacing_impact; return;;
		--help)
			usage spacing_impact; return;;
		*)
			echo "Error in arguments"
			echo $1; usage spacing_impact; return;;
	esac
done

local Errors=0
if [ -z $names ]; then echo "ERROR: -n argument needed or need to be a LIST";Errors+=1;fi
if [ -z $comp_dir ]; then echo "ERROR: -cp argument needed";Errors+=1;fi
if [ -z $out_dir ]; then echo "ERROR: -od argument needed";Errors+=1;fi
if [ -z $genome ]; then echo "ERROR: -g argument needed";Errors+=1;fi
if [ -z $matrix ]; then echo "ERROR: -m argument needed";Errors+=1;fi
if [ -z $spacing ]; then echo "ERROR: -sp argument needed or need to be a LIST";Errors+=1;fi

if [ -z $maxSpace ]; then echo "-maxs argument not used, computing up to 50bp spacing";local maxSpace=50 ;fi
if [ -z $minSpace ]; then echo "-mins argument not used, starting with spacing 0";local minSpace=0 ;fi
if [ -z $offset_left ]; then echo "-maxs argument not used, no left offset used";local offset_left=0 ;fi
if [ -z $offset_right ]; then echo "-maxs argument not used, no right offset used";local offset_right=0 ;fi
if [ -z $colors ]; then echo "-c argument not used, using default palette or need to be a LIST"; local colors=('#40A5C7' '#307C95' '#205364');fi
# TODO check for error in th
if [ $Errors -gt 0 ]; then usage spacing_impact; return 1; fi

mkdir -p -m 774 $out_dir
sed "1d" ${comp_dir}/${names[0]}_${names[1]}/table_${names[0]}_${names[1]}.csv > $out_dir/table_${names[0]}_${names[1]}.bed
local table=$out_dir/table_${names[0]}_${names[1]}.bed

awk -v OFS="\t" '{print $1,$2,$3}' $table > $out_dir/${names[0]}_${names[1]}.bed
local peaks=$out_dir/${names[0]}_${names[1]}.bed

bedtools getfasta -fi $genome -fo $out_dir/${names[0]}_${names[1]}.fa -bed $peaks
local peak_file_fas=$out_dir/${names[0]}_${names[1]}.fa

# python ${scores_prog} -m ${matrix} -f $peak_file_fas -o $out_dir/

# $Python_TFFM $tffm_all_scores -o $out_dir/${names[0]}_${names[1]}_tffm_scores_pos.tsv -pos $out_dir/${names[0]}_${names[1]}.fa -t ${tffm}

python $spacing_mk -o $out_dir/${names[0]}_${names[1]} -smax $maxSpace -smin $minSpace -pos $out_dir/${names[0]}_${names[1]}.fa.scores -th $starting_th -ol $offset_left -or $offset_right -lm 0

awk 'NR==1{save=""; next}NR==2{prechead=$1":"$2":"$3;save=$4$5"\t"$6"\t"$7;next}{head=$1":"$2":"$3;if(head==prechead){save=save"\t"$4$5"\t"$6"\t"$7};if(head!=prechead){print prechead,save;save=$4$5"\t"$6"\t"$7};prechead=head}END{print prechead,save}' $out_dir/${names[0]}_${names[1]}_spacing.tsv > $out_dir/${names[0]}_${names[1]}_spacing.tmp



local spacing=`join_by "," "${spacing[@]}"`
echo $spacing
local files=()
local th1=$starting_th
local continu1=1

while [ $continu1 -eq 1 ];
do
	echo -e "th\t$th1"
	sed 's/[:,/]/\t/g' $out_dir/${names[0]}_${names[1]}_spacing.tmp | awk -v spacing=$spacing -v th1=$th1 -v th2=$th1 -v OFS="\t" '{save=0; split(spacing,list_spacing,","); for(i=4;i<=NF;i+=3){j=i+1;k=i+2;if(($j>=th1 && $k>=th2)||($j>=th2 && $k>=th1)){gsub(/^ER|^IR|^DR|^NA/,"",$i);for(h in list_spacing){if($i == list_spacing[h]){save=1}}}};print $1,$2,$3,save}' | awk -v FS="\t" '{print $4}' > $out_dir/Spacing${th1}.inter

	paste $table $out_dir/Spacing${th1}.inter | sed "1ichr\tbegin\tend\t${names[0]}\t${names[1]}\tSpacing" >  $out_dir/table_${th1}.csv
	files+=("$out_dir/table_${th1}.csv")
	local th1=$(calc $th1+$increment_th)
	local continu1=$(awk -vn1="$th1" -vn2=$ending_th 'BEGIN{print (n1>=n2)?1:0}')
	
done

local files=`join_by "," "${files[@]}"`
local colors=`join_by "," "${colors[@]}"`
$R_36 $impactspacing_mk -f $files -n "${names[0]},${names[1]}" -d $out_dir -c $colors
}



# TODO @Romain, a garder cette version?
compute_rpkmrip_rpkmril_newV(){
	# TODO RBM, corriger cette fonction qui empeche la bonne normalization des données lors du passage RiPvsRiL
	local getCov="yes"; 
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "file containing peaks set to ${2}";shift 2;;
			-bd)
				local bam_dir=("${!2}")
				echo "path to mapping directory set to ${bam_dir[@]}";shift 2;;
			-pd)
				local peaks_dir=("${!2}")
				echo "path to peaks calling directory set to ${peaks_dir[@]}";shift 2;;
			-sn)
				local samples_names=("${!2}")
				echo "list of sample's name set to: ${samples_names[@]}";shift 2;;
			-m)
				local mode=$2
				echo "normalization mode set to 'inPeaks' or 'inLibs' (number of reads retained by MACS2) ${2}";shift 2;;
			-o)
				local out_dir=$2
				echo "out directory set to : ${2}";shift 2;;
			-g)
				local getCov=$2
				echo "optional, choose 'no' if you have already got the reads count (bedtools step): ${2}";shift 2;;
			-h)
				usage compute_rpkmrip_rpkmril; return;;
			--help)
				usage compute_rpkmrip_rpkmril; return;;
			*)
				echo "Error in arguments"
				echo $1; usage compute_rpkmrip_rpkmril; return;;
		esac	
	done
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
	if [ -z $bam_dir ]; then echo "ERROR: -bd argument needed"; Errors+=1; fi
	if [ -z $peaks_dir ]; then echo "ERROR: -pd argument needed"; Errors+=1; fi
	if [ -z $samples_names ]; then echo "ERROR: -sn argument needed or need to be a LIST"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
	if [ -z $mode ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
	if [ -z $getCov ]; then echo "reads count will be computed using bedtools"; fi
	if [ $Errors -gt 0 ]; then usage compute_rpkmrip; return 1; fi	

	mkdir -p -m 774 $out_dir
	#rm $out_dir/tmpFiltTags.txt
	touch $out_dir/tmpTotalTags.txt

	awk -v OFS="\t" 'NR==1 && $2!="begin" && $2!="start"{print $1,$2,$3}NR>1{print $1,$2,$3}' $peaks | sort -k1,1 -k2,2n > $out_dir/tmp_peaks

	printf "pr -mts <(cut -f -3 $out_dir/tmp_peaks) " > $out_dir/tmp_run.sh
	# cut -f -3 $peaks > $out_dir/tmp_peaks # DEPRECATED JL 09/11/2021

	if [[ -f $out_dir/tmpTotalTags.txt ]]; then
		rm $out_dir/tmpTotalTags.txt
		touch $out_dir/tmpTotalTags.txt
	fi
	if [[ -f $out_dir/tmpFiltTags.txt ]]; then
		rm $out_dir/tmpFiltTags.txt
		touch $out_dir/tmpFiltTags.txt
	fi
	## same thing for FRIP file
	if [[ -f ${out_dir}/RIP.txt ]]; then
		rm ${out_dir}/RIP.txt
	fi


			
	for SAMP in  ${samples_names[@]}
	do
		printf "\n\ncomputing reads coverage per peak using bedtools coverage...\n\n"
		echo $SAMP
		printf "<(cut -f 4 ${out_dir}/${SAMP}_filt.cov.bed) " >> $out_dir/tmp_run.sh
		
		for Peaksdir in ${peaks_dir[@]}
		do
			echo $Peaksdir
			if [ -f $Peaksdir/${SAMP}/${SAMP}_stats.txt ]; then
				cut -d " " -f 2 $Peaksdir/$SAMP/${SAMP}_stats.txt | grep -v totalTags >> $out_dir/tmpTotalTags.txt
				cut -d " " -f 3 $Peaksdir/${SAMP}/${SAMP}_stats.txt | grep -v filtTags >> $out_dir/tmpFiltTags.txt
				#break
			fi
		done
		for bams in ${bam_dir[@]}; do
			echo $bams
			
			local bam_file=$(find ${bams} -name "$SAMP.filtered.sorted.dedup.bam" -type f )
			echo "bam_file"
			echo $bam_file
			
			if [ $bam_file != "" ]; then break; fi
		done
		
		if [ "$getCov" == "yes" ]
		then
		bedtools coverage \
			-a $out_dir/tmp_peaks \
			-b $bam_file -sorted -F 1 > ${out_dir}/${SAMP}_filt.cov.bed
		else
			printf "\n\nskipping bedtools coverage because reads coverage has already been computed\n\n"
		fi
		
		inpeaks=$(bedtools sort -i $out_dir/tmp_peaks | bedtools merge -i stdin | bedtools intersect -u -a $bam_file -b stdin -ubam | $PATH_TO_SAMTOOLS/samtools view -c)
		total=$($PATH_TO_SAMTOOLS/samtools view -c $bam_file) # this is the same as what's in tmpTotalTags
		FRIP=$(calc $inpeaks/$total*100)
		echo "${FRIP}% of reads in peaks (${inpeaks}/${total})"
		echo "${inpeaks}" >> ${out_dir}/RIP.txt
			
		#12 aout 2022
		local scale=$(calc 1000000/$inpeaks)
		bedtools genomecov -bga -scale $scale -ibam $bam_file | sed 's/Chr/chr/g' > ${out_dir}/${SAMP}_cpmrip.bdg
		
		#generate bigwig
		cat /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas.fai > ${out_dir}/tmp.chromsize
		$bdg2bwig  ${out_dir}/${SAMP}_cpmrip.bdg ${out_dir}/tmp.chromsize ${out_dir}/${SAMP}_cpmrip.bigwig
		rm ${out_dir}/${SAMP}_cpmrip.bdg
	done

	printf "> ${out_dir}/allReps_RC.txt" >> $out_dir/tmp_run.sh
	bash $out_dir/tmp_run.sh
	rm $out_dir/tmp_run.sh
	#cp ${out_dir}/allReps_RC.txt .
	printf "\n\n\n\nR script to normalize reads count in peaks or in library/mapped running...\n\n"

	/home/prog/R/R-4.4.0/el8/bin/Rscript $rpkmrip_rpkmril ${out_dir}/allReps_RC.txt $mode $out_dir ${samples_names[@]} #TODO replace R-3-5 by R3.6
	# tmpTotalTags is used instead of tmpFiltTags for inLibs normalization in the cmd above

	printf "\n\nend of rpkmrip_rpkmril\n\n"
}

compute_rpkmrip_rpkmril_BK_26nov2024(){ #_BK_19aout2022
	# remise en fonction par JL le 20/10/2022 car la version plus récente créer des pbl
	local getCov="yes"; 
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "file containing peaks set to ${2}";shift 2;;
			-bd)
				local bam_dir=("${!2}")
				echo "path to mapping directory set to ${bam_dir[@]}";shift 2;;
			-pd)
				local peaks_dir=("${!2}")
				echo "path to peaks calling directory set to ${peaks_dir[@]}";shift 2;;
			-sn)
				local samples_names=("${!2}")
				echo "list of sample's name set to: ${samples_names[@]}";shift 2;;
			-m)
				local mode=$2
				echo "normalization mode set to 'inPeaks' or 'inLibs' (number of reads retained by MACS2) ${2}";shift 2;;
			-o)
				local out_dir=$2
				echo "out directory set to : ${2}";shift 2;;
			-g)
				local getCov=$2
				echo "optional, choose 'no' if you have already got the reads count (bedtools step): ${2}";shift 2;;
			-h)
				usage compute_rpkmrip_rpkmril; return;;
			--help)
				usage compute_rpkmrip_rpkmril; return;;
			*)
				echo "Error in arguments"
				echo $1; usage compute_rpkmrip_rpkmril; return;;
		esac	
	done
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
	if [ -z $bam_dir ]; then echo "ERROR: -bd argument needed"; Errors+=1; fi
	if [ -z $peaks_dir ]; then echo "ERROR: -pd argument needed"; Errors+=1; fi
	if [ -z $samples_names ]; then echo "ERROR: -sn argument needed or need to be a LIST"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
	if [ -z $mode ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
	if [ -z $getCov ]; then echo "reads count will be computed using bedtools"; fi
	if [ $Errors -gt 0 ]; then usage compute_rpkmrip; return 1; fi	

	mkdir -p -m 774 $out_dir
	#rm $out_dir/tmpFiltTags.txt
	touch $out_dir/tmpTotalTags.txt

	awk -v OFS="\t" 'NR==1 && $2!="begin" && $2!="start"{print $1,$2,$3}NR>1{print $1,$2,$3}' $peaks | sort -k1,1 -k2,2n > $out_dir/tmp_peaks

	printf "pr -mts <(cut -f -3 $out_dir/tmp_peaks) " > $out_dir/tmp_run.sh
	# cut -f -3 $peaks > $out_dir/tmp_peaks # DEPRECATED JL 09/11/2021

	if [ getCov=="yes" ]; then 
		## if tmpFiltTags file is already present, remove it
		if [[ -f $out_dir/tmpTotalTags.txt ]]; then
			rm $out_dir/tmpTotalTags.txt
			touch $out_dir/tmpTotalTags.txt
		fi
		if [[ -f $out_dir/tmpFiltTags.txt ]]; then
			rm $out_dir/tmpFiltTags.txt
			touch $out_dir/tmpFiltTags.txt
		fi
		## same thing for FRIP file
		if [[ -f ${out_dir}/RIP.txt ]]; then
			rm ${out_dir}/RIP.txt
		fi
	fi


	for SAMP in  ${samples_names[@]}
	do
		printf "\n\ncomputing reads coverage per peak using bedtools coverage...\n\n"
		echo $SAMP
		printf "<(cut -f 4 ${out_dir}/${SAMP}_filt.cov.bed) " >> $out_dir/tmp_run.sh
		
		
		for bams in ${bam_dir[@]}; do
			echo $bams
			
			if find ${bams} -name "$SAMP.filtered.sorted.dedup.bam" -type f -exec false {} +
			then
				local bam_file=$(find ${bams} -name "$SAMP.filtered.sorted.bam" -type f )
			else
				local bam_file=$(find ${bams} -name "$SAMP.filtered.sorted.dedup.bam" -type f )
			fi


			# local bam_file=$(find ${bams} -name "$SAMP.filtered.sorted.bam" -type f )
			echo "bam_file"
			echo $bam_file
			
			if [ $bam_file != "" ]; then 
				break
			fi
			
			
		done
		#  this is to remove the control bam # NOTE NEW by JL 02/12/2021
		
		# # check that the bam file listed is a file and not a link, and
		# # remove the link from the array if present
		# for bam in ${bam_file[@]}; do
		# 	if [[ ! -L $bam && -f $bam ]];then
		# 		echo "bam file is the correct one: $bam"
		# 	else
		# 		echo "bam_file is a link, we do not consider it"
		# 		bam_file=("${bam_file[@]/$bam}")
				
		# 	fi
		# done
		
		
		# DEPRECATED by JL 02/12/2021
		# local tmp241=$bam_dir/$SAMP/*.filtered.sorted.nodup.bam
		# echo ${tmp241}
		# local bam_file=$(echo $tmp241 | sed -e 's/ .*$//g') #this is to remove the control bam
		echo "hello"
		echo ${bam_file}
		if [ "$getCov" == "yes" ]
		then
		bedtools coverage \
			-a $out_dir/tmp_peaks \
			-b $bam_file -F 0.5 > ${out_dir}/${SAMP}_filt.cov.bed
			inpeaks=$(bedtools sort -i $out_dir/tmp_peaks | bedtools merge -i stdin | bedtools intersect -u -a $bam_file -b stdin -ubam | $PATH_TO_SAMTOOLS/samtools view -c)
			total=$($PATH_TO_SAMTOOLS/samtools view -c $bam_file)
			FRIP=$(calc $inpeaks/$total*100)
			echo "${FRIP}% of reads in peaks (${inpeaks}/${total})"
			echo "${inpeaks}" >> ${out_dir}/RIP.txt
			
			#12 aout 2022
			local scale=$(calc 1000000/$inpeaks)
			bedtools genomecov -bga -scale $scale -ibam $bam_file | sed 's/Chr/chr/g' > ${out_dir}/${SAMP}_cpmrip.bdg
			
			#format filtered tags (for in libs normalization)
			
			for Peaksdir in ${peaks_dir[@]}
			do
				if [ -f $peaks_dir/$SAMP/${SAMP}_stats.txt ]; then
					cut -d " " -f 2 $peaks_dir/$SAMP/${SAMP}_stats.txt | grep -v totalTags >> $out_dir/tmpTotalTags.txt
					cut -d " " -f 3 $peaks_dir/$SAMP/${SAMP}_stats.txt | grep -v filtTags >> $out_dir/tmpFiltTags.txt
					break
				fi
			done
		else
			printf "\n\nskipping bedtools coverage because reads coverage has already been computed\n\n"
		fi
		
	done

	printf "> ${out_dir}/allReps_RC.txt" >> $out_dir/tmp_run.sh
	bash $out_dir/tmp_run.sh
	rm $out_dir/tmp_run.sh
	#cp ${out_dir}/allReps_RC.txt .
	printf "\n\n\n\nR script to normalize reads count in peaks or in library/mapped running...\n\n"
	/home/prog/R/R-4.4.0/el8/bin/Rscript $rpkmrip_rpkmril ${out_dir}/allReps_RC.txt $mode $out_dir ${samples_names[@]} #TODO replace R-3-5 by R3.6
	printf "\n\nend of rpkmrip_rpkmril\n\n"
}

compute_rpkmrip_rpkmril(){ 
	# remise en fonction par JL le 20/10/2022 car la version plus récente créer des pbl
	local getCov="yes"; 
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "file containing peaks set to ${2}";shift 2;;
			-bd)
				local bam_dir=("${!2}")
				echo "path to mapping directory set to ${bam_dir[@]}";shift 2;;
			-pd)
				local peaks_dir=("${!2}")
				echo "path to peaks calling directory set to ${peaks_dir[@]}";shift 2;;
			-sn)
				local samples_names=("${!2}")
				echo "list of sample's name set to: ${samples_names[@]}";shift 2;;
			-m)
				local mode=$2
				echo "normalization mode set to 'inPeaks' or 'inLibs' (number of reads retained by MACS2) ${2}";shift 2;;
			-o)
				local out_dir=$2
				echo "out directory set to : ${2}";shift 2;;
			-g)
				local getCov=$2
				echo "optional, choose 'no' if you have already got the reads count (bedtools step): ${2}";shift 2;;
			-h)
				usage compute_rpkmrip_rpkmril; return;;
			--help)
				usage compute_rpkmrip_rpkmril; return;;
			*)
				echo "Error in arguments"
				echo $1; usage compute_rpkmrip_rpkmril; return;;
		esac	
	done
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument needed"; Errors+=1; fi
	if [ -z $bam_dir ]; then echo "ERROR: -bd argument needed"; Errors+=1; fi
	if [ -z $peaks_dir ]; then echo "ERROR: -pd argument needed"; Errors+=1; fi
	if [ -z $samples_names ]; then echo "ERROR: -sn argument needed or need to be a LIST"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -o argument needed"; Errors+=1; fi
	if [ -z $mode ]; then echo "ERROR: -m argument needed"; Errors+=1; fi
	if [ -z $getCov ]; then echo "reads count will be computed using bedtools"; fi
	if [ $Errors -gt 0 ]; then usage compute_rpkmrip; return 1; fi	

	mkdir -p -m 774 $out_dir
	#rm $out_dir/tmpFiltTags.txt
	touch $out_dir/tmpTotalTags.txt

	awk -v OFS="\t" 'NR==1 && $2!="begin" && $2!="start"{print $1,$2,$3}NR>1{print $1,$2,$3}' $peaks | sort -k1,1 -k2,2n > $out_dir/tmp_peaks

	printf "pr -mts <(cut -f -3 $out_dir/tmp_peaks) " > $out_dir/tmp_run.sh
	# cut -f -3 $peaks > $out_dir/tmp_peaks # DEPRECATED JL 09/11/2021

	if [ getCov=="yes" ]; then 
		## if tmpFiltTags file is already present, remove it
		if [[ -f $out_dir/tmpTotalTags.txt ]]; then
			rm $out_dir/tmpTotalTags.txt
			touch $out_dir/tmpTotalTags.txt
		fi
		if [[ -f $out_dir/tmpFiltTags.txt ]]; then
			rm $out_dir/tmpFiltTags.txt
			touch $out_dir/tmpFiltTags.txt
		fi
		## same thing for FRIP file
		if [[ -f ${out_dir}/RIP.txt ]]; then
			rm ${out_dir}/RIP.txt
		fi
	fi


	for SAMP in  ${samples_names[@]}
	do
		printf "\n\ncomputing reads coverage per peak using bedtools coverage...\n\n"
		echo $SAMP
		printf "<(cut -f 4 ${out_dir}/${SAMP}_filt.cov.bed) " >> $out_dir/tmp_run.sh
		
		bam_file="NA"
		for bams in ${bam_dir[@]}; do
			echo $bams
			
			if [[ -f "${bams}/$SAMP/$SAMP.filtered.sorted.dedup.bam" ]]; then
				local bam_file=${bams}/$SAMP/$SAMP.filtered.sorted.dedup.bam
				break
			elif [[ -f "${bams}/$SAMP/$SAMP.filtered.sorted.bam" ]]; then
				local bam_file=${bams}/$SAMP/$SAMP.filtered.sorted.bam
				break
			fi

			echo "bam_file"
			echo $bam_file
			if [[ -f "$bam_file" ]]; then
				break
			fi
			
		done
		echo ${bam_file}

		if [ "$getCov" == "yes" ]
		then
		# bedtools coverage \
			# -a $out_dir/tmp_peaks \
			# -b $bam_file -F 0.5 > ${out_dir}/${SAMP}_filt.cov.bed
			echo "bed cov"
			bedtools coverage \
			-a $out_dir/tmp_peaks \
			-b $bam_file > ${out_dir}/${SAMP}_filt.cov.bed
			inpeaks=$(bedtools sort -i $out_dir/tmp_peaks | bedtools merge -i stdin | bedtools intersect -u -a $bam_file -b stdin -ubam | $PATH_TO_SAMTOOLS/samtools view -c)
			total=$($PATH_TO_SAMTOOLS/samtools view -c $bam_file)
			FRIP=$(calc $inpeaks/$total*100)
			echo "${FRIP}% of reads in peaks (${inpeaks}/${total})"
			echo "${inpeaks}" >> ${out_dir}/RIP.txt
			#12 aout 2022
			# local scale=$(calc 1000000/$inpeaks)
			# local scale=1
			# bedtools genomecov -bga -scale $scale -ibam $bam_file | sed 's/Chr/chr/g' > ${out_dir}/${SAMP}_cpmrip.bdg
			
			#format filtered tags (for in libs normalization)
			local found=0
			for Peaksdir in ${peaks_dir[@]}
			do
				if [ -f $peaks_dir/$SAMP/${SAMP}_stats.txt ]; then
					cut -d " " -f 2 $peaks_dir/$SAMP/${SAMP}_stats.txt | grep -v totalTags >> $out_dir/tmpTotalTags.txt
					cut -d " " -f 3 $peaks_dir/$SAMP/${SAMP}_stats.txt | grep -v filtTags >> $out_dir/tmpFiltTags.txt
					local found=1
					break
				fi
			done
			echo $found $bam_file
			if [ $found -eq 0 ] && [[ $bam_file == *"control"* ]]; then
				echo $(($total)) >> $out_dir/tmpTotalTags.txt
				echo $(($total)) >>  $out_dir/tmpFiltTags.txt
			fi
		else
			printf "\n\nskipping bedtools coverage because reads coverage has already been computed\n\n"
		fi
	done
	
	printf "> ${out_dir}/allReps_RC.txt" >> $out_dir/tmp_run.sh
	bash $out_dir/tmp_run.sh
	rm $out_dir/tmp_run.sh
	#cp ${out_dir}/allReps_RC.txt .
	printf "\n\n\n\nR script to normalize reads count in peaks or in library/mapped running...\n\n"
		echo '~###################################'

	/home/prog/R/R-4.4.0/el8/bin/Rscript $rpkmrip_rpkmril ${out_dir}/allReps_RC.txt $mode $out_dir ${samples_names[@]} #TODO replace R-3-5 by R3.6
	printf "\n\nend of rpkmrip_rpkmril\n\n"
}

initial_comparison(){
	# initial_comparison -n1 -n2 -od -id -bd -p1 -p2
	local get_bedtools_cov="yes"; local filterCov=0; local filterHeight=0; local bamdir2="NA"; local data2="NA"; local peaksext1="NA"; local peaksext2="NA"
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-n1)
				local name1=$2
				echo "name of directory for dataset 1 is: ${2}";shift 2;;
			-n2)
				local name2=$2
				echo "name of directory for dataset 2 is: ${2}";shift 2;;
			-od)
				local result=$2
				echo "output directory set to: ${2}";shift 2;;
			-id)
				local data=$2
				echo "general data directory (i.e. peakcalling directory) set to: ${2}";shift 2;;
			-id2)
				local data2=$2
				echo "general data directory (i.e. peakcalling directory) for second dataset set to: ${2}";shift 2;;
			-f)
				local filterCov=$2
				echo "coverage must be at least >${2} in both sample to be considered as peaks";shift 2;;
			-he)
				local filterHeight=$2
				echo "height must be at least >${2} in both sample to be considered as peaks";shift 2;;
			-bd)
				local bamdir=$2
				echo "general bam directory set to: ${2}";shift 2;;
			-p1)
				local peaksext1=$2
				echo "peaks to use set to: ${2}";shift 2;;
			-p2)
				local peaksext2=$2
				echo "peaks to use set to: ${2}";shift 2;;
			-bd2)
				local bamdir2=$2
				echo "general bam directory for second dataset set to: ${2}";shift 2;;
			-gcov)
				local get_bedtools_cov=$2
				echo "'yes' or 'no' to compute reads count from bam at each peak: ${2}";shift 2;;
			-rep1)
				local list_rep1=("${!2}")
				echo "list of replicates names for dataset 1 in bam directory set to: ${list_rep1[@]}";shift 2;;
			-rep2)
				local list_rep2=("${!2}")
				echo "list of replicates names for dataset 2 in bam directory set to: ${list_rep2[@]}";shift 2;;
			-h)
				usage initial_comparison; return;;
			--help)
				usage initial_comparison; return;;
			*)
				echo "Error in arguments"
				echo $1; usage initial_comparison; return;;
		esac
	done
	local Errors=0
	if [ -z $name1 ]; then echo "ERROR: -n1 argument needed"; Errors+=1; fi
	if [ -z $name2 ]; then echo "ERROR: -n2 argument needed"; Errors+=1; fi
	if [ -z $result ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ -z $data ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage initial_comparison; return 1; fi
	if [ $data2 == "NA" ]; then
		echo "id2 not used, assuming same Peakcalling directory for both datasets"
		local data2=$data
	fi
	if [ $bamdir2 == "NA" ]; then
		echo "bd2 not used, assuming same bam directory for both datasets"
		local bamdir2=$bamdir
	fi


	# export TMPDIR=/nobackup
	# local do_plot_quant=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/Hist_cov_gen.r
	# local merge_peaks=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/merge_all_peaks.py
	# local compute_coverage=/home/312.6-Flo_Re/312.6.1-Commun/scripts/DAP_global_analysis_p3.7/compute_coverage.py

	local out_dir=$result/${name1}_${name2}
	mkdir -p -m 774 $out_dir
	local cov1=$data/$name1/${name1}_cov.bdg
	local cov2=$data2/$name2/${name2}_cov.bdg
	local remove_list=()
	###### few lines added for selecting peaks on their height:
	if [ $peaksext1 == "NA" ]; then
		if [ $filterHeight -eq 0 ];then
			local peaks1=$data/$name1/${name1}_narrow.bed
		else
			awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data/$name1/${name1}_max.bed > $data/$name1/${name1}_narrow.heightFiltered.bed
			local peaks1=$data/$name1/${name1}_narrow.heightFiltered.bed
		fi
	else
		#TODO allow filter height here ?
		local peaks1=$peaksext1
	fi
	if [ $peaksext2 == "NA" ]; then
		if [ $filterHeight -eq 0 ];then
			local peaks2=$data2/$name2/${name2}_narrow.bed
		else
			awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data2/$name2/${name2}_max.bed > $data2/$name2/${name2}_narrow.heightFiltered.bed
			local peaks2=$data2/$name2/${name2}_narrow.heightFiltered.bed
		fi
	else
		#TODO allow filter height here ?
		local peaks2=$peaksext2
	fi
 
		# retired 20/03/2024 by JL to autorise the use of external peaks for 1 or both dataset
		# if [ $filterHeight -eq 0 ];then
		# 	local peaks1=$data/$name1/${name1}_narrow.bed
		# 	local peaks2=$data2/$name2/${name2}_narrow.bed
		# else
		# 		# awk -v filt=$filterHeight -v OFS="\t" '$4>filt && $5>filt && $6>filt' $data/$name1/${name1}_max.bed > $data/$name1/${name1}_narrow.heightFiltered.bed
		# 		awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data/$name1/${name1}_max.bed > $data/$name1/${name1}_narrow.heightFiltered.bed
		# 		local peaks1=$data/$name1/${name1}_narrow.heightFiltered.bed
		# 		# awk -v filt=$filterHeight -v OFS="\t" '$4>filt && $5>filt && $6>filt' $data/$name2/${name2}_max.bed > 	$data/$name2/${name2}_narrow.heightFiltered.bed
		# 		awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data2/$name2/${name2}_max.bed > $data2/$name2/${name2}_narrow.heightFiltered.bed
		# 		local peaks2=$data2/$name2/${name2}_narrow.heightFiltered.bed
		# fi


    awk -v OFS="\t" -v name=$name1 '{print $1,$2,$3,name}' $peaks1  > $out_dir/${name1}_peaks.bed
    awk -v OFS="\t" -v name=$name2 '{print $1,$2,$3,name}' $peaks2  > $out_dir/${name2}_peaks.bed
    
    cat $out_dir/${name1}_peaks.bed $out_dir/${name2}_peaks.bed | sort -k1,1 -k2,2n > $out_dir/${name1}_${name2}_peaks.bed
    python3 $merge_peaks -f1 $name1 -f2 $name2 -o $result
    echo "merged"
    local peak_file=$out_dir/${name1}_${name2}_peaks_processed.bed
    cat $out_dir/${name1}_${name2}_peaks_merged.bed $out_dir/${name1}_${name2}_peaks_uniques.bed | sort -k1,1 -k2,2n | uniq > $peak_file
    echo "sorted"
    bedtools intersect -a $peak_file -b $cov1 -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > $out_dir/$name1.inter &
    bedtools intersect -a $peak_file -b $cov2 -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > $out_dir/$name2.inter &
    wait
    echo "intersected"
    python3 $compute_coverage -i $out_dir/$name1.inter -m
    python3 $compute_coverage -i $out_dir/$name2.inter -m
	remove_list+=("$out_dir/$name1.inter" "$out_dir/$name2.inter")
    echo "cov computed"
    awk  '{print (1000*$5)/($3-$2)}' $out_dir/$name1.inter.cov  > $out_dir/$name1.inter.cov.normed # adjust cov by len peak
    awk  '{print (1000*$5)/($3-$2)}' $out_dir/$name2.inter.cov  > $out_dir/$name2.inter.cov.normed
	# adjust cov by len peak
	# awk  '($3-$2)!=0{print (1000*$5)/($3-$2);next}($3-$2)==0{print (1000*$5)/1}' $out_dir/$name1.inter.cov  > $out_dir/$name1.inter.cov.normed # JL 19/01/2023 added control for 
    # awk  '($3-$2)!=0{print (1000*$5)/($3-$2);next}($3-$2)==0{print (1000*$5)/1}' $out_dir/$name2.inter.cov  > $out_dir/$name2.inter.cov.normed
    echo "norm adjusted"
	if [ $filterCov -eq 0 ];then
		paste $peak_file $out_dir/$name1.inter.cov.normed $out_dir/$name2.inter.cov.normed  |  sed "1ichr\tbegin\tend\tname\t$name1\t$name2" >  $out_dir/table_${name1}_${name2}.csv
	else
		paste $peak_file $out_dir/$name1.inter.cov.normed $out_dir/$name2.inter.cov.normed | awk -v filter=$filterCov -v OFS="\t" '$5>=filter || $6>=filter {print $0}' |  sed "1ichr\tbegin\tend\tname\t$name1\t$name2" >  $out_dir/table_${name1}_${name2}.csv
	fi
	
	if [ ! -z $bamdir ] && [ -z $list_rep1 ]; then # Bamdir defined but not list_rep
		local list_rep1=()
		local list_rep2=()
		for rep in {1..10}; do ## LT
			echo $rep ## LT
			# echo $bamdir/${name1}rep${rep} ## LT
			if [[ -d $bamdir/${name1}rep${rep} ]]; then
				list_rep1+=("${name1}rep${rep}")
			fi
			if [[ -d $bamdir/${name1}_rep${rep} ]]; then
				list_rep1+=("${name1}_rep${rep}")
			fi

			# echo $bamdir/${name2}rep${rep} ## LT
			if [[ -d $bamdir2/${name2}rep${rep} ]]; then
				list_rep2+=("${name2}rep${rep}")
			fi
			if [[ -d $bamdir2/${name2}_rep${rep} ]]; then
				list_rep2+=("${name2}_rep${rep}")
			fi
			if [[ ! -d $bamdir/${name1}rep${rep} ]] && [[ ! -d $bamdir2/${name2}rep${rep} ]]; then
				break 
				# continue ## LT
			fi
		done
	# exit 0
	
	fi
	if [ -z $bamdir ] && [ -z $list_rep1 ]; then # Bamdir & list_rep not defined
		/home/prog/R/R-4.4.0/el8/bin/Rscript $do_plot_quant $result $name1 $name2 $out_dir/table_${name1}_${name2}.csv RiL #TODO replace R-3-5 by R3.6
	else
		echo ${list_rep1[@]}
		echo ${list_rep2[@]}
		local list_all=("${list_rep1[@]}" "${list_rep2[@]}")
		echo ${list_all[@]}
		local list_bamdir=("$bamdir" "$bamdir2")
		local list_data=("$data" "$data2")
		compute_rpkmrip_rpkmril -p $out_dir/table_${name1}_${name2}.csv -bd list_bamdir[@] -pd list_data[@] -sn list_all[@] -m "inLibs" -o $out_dir/ -g $get_bedtools_cov
		last_col_n1="$((${#list_rep1[@]}-1+4))"
		echo $last_col_n1
		first_col_n2="$(($last_col_n1+1))"
		echo $first_col_n2
		awk -v OFS="\t" -v c=$last_col_n1 'NR==1 && $2!="begin" && $2!="start"{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}NR>1{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name1}_RPKMril.txt
		awk -v OFS="\t" -v c=$first_col_n2 'NR==1 && $2!="begin" && $2!="start"{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}NR>1{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminLibs.txt > $out_dir/${name2}_RPKMril.txt
		
		compute_rpkmrip_rpkmril -p $out_dir/table_${name1}_${name2}.csv -bd list_bamdir[@] -pd list_data[@] -sn list_all[@] -m "inPeaks" -o $out_dir/ -g "no"
		awk -v OFS="\t" -v c=$last_col_n1 'NR==1 && $2!="begin" && $2!="start"{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}NR>1{moy=0;nb=0;for(i=4;i<=c;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminPeaks.txt > $out_dir/${name1}_RPKMrip.txt
		awk -v OFS="\t" -v c=$first_col_n2 'NR==1 && $2!="begin" && $2!="start"{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}NR>1{moy=0;nb=0;for(i=c;i<=NF;i++){nb++;moy+=$i};print $1,$2,$3, (moy)/nb}' $out_dir/peaks_perSample_rpkminPeaks.txt > $out_dir/${name2}_RPKMrip.txt

		if [ -f $out_dir/${name1}_RPKMrip.txt ] && [ -f $out_dir/${name2}_RPKMrip.txt ] && [ -f $out_dir/${name1}_RPKMril.txt ] && [ -f $out_dir/${name2}_RPKMril.txt ]; then
			echo "name1 : ${name1}"
			echo "name2 : ${name2}"
			if [ ${#list_rep1[@]} -gt 1 ] && [ ${#list_rep2[@]} -gt 1 ]; then

				/home/prog/R/R-4.4.0/el8/bin/Rscript /home/312.6-Flo_Re/312.6.1-Commun/ARF-anr/DAP_052022/statsedgeR.R $out_dir/peaks_perSample_rpkminPeaks.txt $out_dir ${#list_rep1[@]} ${#list_rep2[@]} 
				paste <(awk 'NR==1 && $2!="begin"{print $0}NR>1{print $0}' $out_dir/table_${name1}_${name2}.csv) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name1}_RPKMrip.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name2}_RPKMrip.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name1}_RPKMril.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name2}_RPKMril.txt) <(awk -v OFS="\t" 'NR!=1 {print $2,$6}' $out_dir/results_edgeR.tsv) | sed "1ichr\tbegin\tend\tname\t${name1}_RiLgenomecov\t${name2}_RiLgenomecov\t${name1}_RiP\t${name2}_RiP\t${name1}_RiL\t${name2}_RiL\tlogFC\tFDR" > $out_dir/table_${name1}_${name2}_RiL_RiP.tsv
			else
				echo "so im not in the loop"
				paste <(awk 'NR==1 && $2!="begin"{print $0}NR>1{print $0}' $out_dir/table_${name1}_${name2}.csv) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name1}_RPKMrip.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name2}_RPKMrip.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name1}_RPKMril.txt) <(awk 'NR==1 && $2!="start"{print $4}NR>1{print $4}' $out_dir/${name2}_RPKMril.txt) | sed "1ichr\tbegin\tend\tname\t${name1}_RiLgenomecov\t${name2}_RiLgenomecov\t${name1}_RiP\t${name2}_RiP\t${name1}_RiL\t${name2}_RiL" > $out_dir/table_${name1}_${name2}_RiL_RiP.tsv
			fi
		fi
		echo "enter R script"
		/home/prog/R/R-4.4.0/el8/bin/Rscript $do_plot_quantRIP $result $name1 $name2 $out_dir/table_${name1}_${name2}_RiL_RiP.tsv RiP #TODO replace R-3-5 by R3.6
		echo "after R script"
		#Rscript $do_plot_quant $result $name1 $name2 $out_dir/table_${name1}_${name2}_RiL_RiP.tsv RiP
		# /home/312.6-Flo_Re/312.6.1-Commun/Programs/Anaconda3/envs/py38_ml/bin/python $Prcurves -m1 ${name1} -m2 ${name2} -fp $out_dir/table_${name1}_${name2}_RiL_RiP.tsv -od $result
	
	fi
	rm $remove_list
	
	# makes the function failed AJ 05/2025 (files do not exist?)
	# if [ ! -z $bamdir ] || [ ! -z $list_rep1 ]; then
	# 	rm $out_dir/*.bdg
	# fi
	
}

main_peakcalling_BK_26nov2024(){
	echo "=== main_peakcalling function started ==="
	local top=0; local mspcpc=100; local phs=200; local Genome_length=120000000; local seedrandom=168159; local additionnal_args=""; local keepDup=1;  local threads=8; local redo_analysis="false"; local input_dir=("NA"); local sizeFile="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.chromsize"; local blacklist="/home/312.6-Flo_Re/312.6.1-Commun/data/A_thaliana_phytozome_v12/Greenscreen_19012023_merged.bed"
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-id)
				local in_dir=("${!2}")
				echo "data directory set to: ${in_dir[@]}";shift 2;;
			-cd)
				local input_dir=("${!2}")
				echo "control directory set to: ${input_dir[@]}";shift 2;;
			-od)
				local out_dir=$2
				echo "output directory set to: ${2}";shift 2;;
			-nc)
				local name_cons=$2
				echo "name for consensus directory set to: ${2}";shift 2;;
			-g)
				local Genome_length=${2%.*}
				echo "genome length (mappable) set to: ${2}";shift 2;;
			-add)
				local additionnal_args="$2" # TODO fix this argument
				echo "additionnal argument set to MACS2: ${2}";shift 2;;
				# CONCERNING previous TODO:
				# it seems this way of giving string argument doesn't really works well
			-top)
				local top=$2
				echo "maximum number of peaks set to: ${2}";shift 2;;
			-ps)
				local phs=$(calc $2/2.0) # $2 is the peak size! but we take only half, hence peak half size...
				echo "peak size: ${2}";shift 2;;
			-size)
				local sizeFile=$2; 
				echo "size file for genome is: $2"; shift 2;;
			-bl)
				local blacklist=$2 ; shift 2 ;;
			-mspc)
				local mspcpc=$2
				echo "percentage for mspc set to: $2"; shift 2;;
			-s)
				local seedrandom=$2
				echo "random seed set to: ${2}"; shift 2;;
			-t)
				local threads=$2; shift 2;;
			--keep-dup)
				local keepDup=$2; shift 2;;
			-r)
				local redo_analysis="true"; shift 1;;
			-h)
				usage main_peakcalling; return;;
			--help)
				usage main_peakcalling; return;;
			*)
				echo "Error in arguments"
				echo $1"\t"$2; usage main_peakcalling; return;;
		esac
	done
	local Errors=0
	if [ -z $in_dir ]; then echo "ERROR: -id argument needed or need to be a LIST"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ -z $input_dir ]; then echo "ERROR: -cd argument needed or need to be a LIST"; Errors+=1; fi
	if [ -z $name_cons ]; then echo "ERROR: -nc argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage main_peakcalling; return 1; fi
	local list_peaks_mspc=()
	local list_bdg=()
	local list_rm=()
	if [[ $mspcpc != *"%"* ]]; then local mspcpc=$mspcpc"%" ; fi
	echo -e "data_dir:$in_dir\ncontrol_dir:$input_dir\noutput_dir:$out_dir\nname_consensus_dir:$name_cons\nGenome_length:$Genome_length\nadditionnal_args:$additionnal_args\nmax_nb_peaks:$top\npeaks_half_size:$phs\nmspc_percentage:$mspcpc\nseed_for_random:$seedrandom\n--keep-dup:$keepDup\n" > $out_dir/${name_cons}/log_params.txt # log for parameters_used
	tmpsave=$LC_COLLATE;LC_COLLATE=C
	for dir in ${in_dir[@]};
	do
		# link bamfile of the controls into the directory of the sample (usefull for later use in MACS2)
		bam_ech=$(find $dir -name "*.filtered.sorted.bam")
		local i=1
		if [ ${input_dir[0]} != "NA" ]; then
			for ctrldir in ${input_dir[@]};
			do
				for bam in $(find $ctrldir -name "*.filtered.sorted.bam")
				do
					if [ ! -f $dir/control_$i.filtered.sorted.bam ] || [[ $bam -nt $dir/control_$i.filtered.sorted.bam ]]; then
						ln -s $bam $dir/control_$i.filtered.sorted.bam
					fi
					local i=$(($i+1))
				done
			done
		fi
		local name=${dir##*/}
		# expected outputs of function peakcalling_MACS2 per replicate
		list_peaks_mspc+=("$out_dir/$name/${name}_peaks.bed") 
		list_bdg+=("$out_dir/$name/${name}_cov.bdg")
		peakcalling_MACS2 -id $dir -od $out_dir -g $Genome_length -s $seedrandom -phs $phs -kd $keepDup -t $threads -size $sizeFile -add "$additionnal_args" -r $redo_analysis -bl $blacklist
		for elt in ${bam_ech[@]};
		do
			if [ ! -f $out_dir/$name/${name}_stats.txt ] || [[ $out_dir/$name/${name}_peaks.bed -nt $out_dir/$name/${name}_stats.txt ]]; then # checking if this step is already done to not do it again
				if [[ $elt != *"control"* ]]; then # filter out all bam that contain control in his name.
					# computing Fresquency of Reads in Peaks (FRiP)
					total=$($PATH_TO_SAMTOOLS/samtools view -c $elt)
					inpeak=$(bedtools sort -i $out_dir/$name/${name}_peaks.bed | bedtools merge -i stdin | bedtools intersect -u -a $elt -b stdin -ubam | $PATH_TO_SAMTOOLS/samtools view -c)
					FRIP=$(calc $inpeak/$total*100)
					echo "$name $FRIP% $inpeak / $total" > $out_dir/$name/${name}_FreqReadInPeak.txt
					R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${elt%.filtered.sorted.bam}.minimal.stats)
					if [ $R2_value == "NA" ]; then # number of tags/fragsment, way of rzetrienving this info is PE/SE dependant
						totTags=$(grep "total tags in treatment" $out_dir/$name/${name}.log | cut -d " " -f 15)
						if [ $keepDup == "all" ]; then
							filtTags=$totTags
						else
							filtTags=$(grep "tags after filtering in treatment" $out_dir/$name/${name}.log | cut -d " " -f 16)
						fi
					fi
					if [ $R2_value != "NA" ]; then
						totTags=$(grep "total fragments in treatment" $out_dir/$name/${name}.log | cut -d " " -f 15)
						if [ $keepDup == "all" ]; then
							filtTags=$totTags
						else
							filtTags=$(grep "fragments after filtering in treatment" $out_dir/$name/${name}.log | cut -d " " -f 16)
						fi
					fi
					filtPeaks=$(wc -l $out_dir/$name/${name}_filtered.narrowPeak | cut -d " " -f 1)
					echo "Sample totalTags filtTags filtPeaks FRIP" > $out_dir/$name/${name}_stats.txt
					echo $name $totTags $filtTags $filtPeaks $FRIP >> $out_dir/$name/${name}_stats.txt
				fi
			fi
		done
	done
	# removing directory if analysis have already been done (starting fresh)
	if [[ -d $out_dir/$name_cons ]]; then rm -Rf $out_dir/$name_cons; fi
	mkdir -p -m 774 $out_dir/$name_cons/temp
	if [ "${#in_dir[@]}" -gt 1 ]; then # if more than one replicates
		# if [[ $(uname -a) =~ "el7" ]]; then # in our sys, el7 is needed to run mspc
	# 			echo "${mspc_mk} -i ${list_peaks_mspc[@]} -r Tec -w 1e-4 -s 1e-8 -c ${mspcpc} -o $out_dir/$name_cons -d 1" #• DEBUG
			# p_val=$(awk 'BEGIN {print 10**-30}') #TODO implement this for MSPC
			${mspc_mk} -i ${list_peaks_mspc[@]} -r Tec -w 1e-4 -s 1e-8 -c ${mspcpc} -o $out_dir/$name_cons -d 1
		# else
		# 	echo "Error SL7 node is needed for mspc" && return 1
		# fi
		# filtering max number of peaks or not, depending on user choice. Filter is done on score
		if [ $top -eq 0 ]; then  
			sed '1d' $out_dir/$name_cons/ConsensusPeaks.bed | awk -v OFS="\t" '{print $1,$2,$3}' | sort -k1,1 -k2,2n > $out_dir/$name_cons/${name_cons}.bed
		else
			sed '1d' $out_dir/$name_cons/ConsensusPeaks.bed | sort -k5,5nr | awk -v OFS="\t" '{print $1,$2,$3}' > $out_dir/$name_cons/${name_cons}_comp.bed
			head -${top} $out_dir/$name_cons/${name_cons}_comp.bed | sort -k1,1 -k2,2n > $out_dir/$name_cons/${name_cons}.bed
		fi
		local i=0
		local files=()
		# prep files for next big filter, here we compute positions for the maximum of each peaks, for each replicate
		for dir in ${in_dir[@]};
		do
			local name=${dir##*/}
			if [ $i -eq 0 ]; then
			# why $5+$13 ? $5 start=filtered narropeak & $13 position of  maximum
				bedtools intersect -a $out_dir/$name_cons/$name_cons.bed -b $out_dir/$name/${name}_filtered.narrowPeak -loj | awk -v OFS="\t" '{print $1,$2,$3,$5+$13}' | awk -v OFS="\t" -v chr="" -v start="" -v stop="" -v save="" 'start!=$2 && stop !=$3{if(save!=""){print save};save=$0;chr=$1;start=$2;stop=$3;next} chr==$1 && start==$2 && stop ==$3 {save=save" "$4}END{print save}' > $out_dir/$name_cons/temp/tmp_peaks_$i.bed 
			else
				bedtools intersect -a $out_dir/$name_cons/$name_cons.bed -b $out_dir/$name/${name}_filtered.narrowPeak -loj | awk -v OFS="\t" '{print $1,$2,$3,$5+$13}' | awk -v OFS="\t" -v chr="" -v start="" -v stop="" -v save="" 'start!=$2 && stop !=$3{if(save!=""){print save};save=$0;chr=$1;start=$2;stop=$3;next} chr==$1 && start==$2 && stop ==$3 {save=save" "$4}END{print save}' | awk -v OFS=" " '{$1=$2=$3="";print $0}' | sed 's/   //' > $out_dir/$name_cons/temp/tmp_peaks_$i.bed
			fi
			files+=("$out_dir/$name_cons/temp/tmp_peaks_$i.bed")
			local i=$(($i+1))
		done
		rm $out_dir/$name_cons/$name_cons.bed # cleaning unused misleading file
		# the maximum of each peaks is compared among each replicates. if maximum are close enough (peak half size) they are considered as same maximum, otherwize multiple peaks are created to insure that 1 maximum => 1 peaks rule is respected.
		paste ${files[@]} > $out_dir/$name_cons/temp/${name_cons}_narrow.bed
		thpc=$(calc ${#files[@]}*${mspcpc//%}/100 )
		echo "$thpc,$phs"
		awk -v SIZE=$phs -v th=${thpc} -v FS="[ \t]" -v OFS="\t" 'function abs(v) {return v < 0 ? -v : v} {nbmax=0; 
		for(i=4;i<=NF;i++) { 
			if($i!= -1){ 
				if(nbmax!=0){
					ok=1;
					for(lmax in max){
						localmax=max[lmax];
						if(abs(localmax-$i)<=SIZE){ 
							moy[localmax]+=$i; nb[localmax]+=1; ok=0; break 
						}
					} 
					if(ok==1){ 
						max[i-3]=$i; moy[$i]=$i;nb[$i]=1
					} 
				} 
				if(nbmax==0){
					max[i-3]=$i; moy[$i]=$i;nb[$i]=1;nbmax=1
				} 
			}
		}; 
		for(kmax in max){
			localmax=max[kmax]; 
			if(nb[localmax]>=th){
				print($1, int(moy[localmax]/nb[localmax])-SIZE, int(moy[localmax]/nb[localmax])+SIZE)
			}; 
			delete max[kmax];delete nb[localmax]; delete moy[localmax]
		}
		}' $out_dir/$name_cons/temp/${name_cons}_narrow.bed | sed 's/\t\+/\t/g;s/^\t//' > $out_dir/$name_cons/${name_cons}_narrow.bed
		# finally the mean coverage values of replicates is used as consensus coverage values
		echo "computing mean coverage for consensus coverage. This might take a while...."
		$PATH_TO_MACS/macs2 cmbreps -i ${list_bdg[@]} -m mean -o $out_dir/$name_cons/${name_cons}_unsorted_cov.bdg
		
	else
		# no step of replicates "merging" is necessary when you have only one replicate
		local name=${in_dir[0]##*/}
		cp $out_dir/$name/${name}_narrow.bed $out_dir/$name_cons/${name_cons}_narrow.bed
		cp ${list_bdg[0]} $out_dir/$name_cons/${name_cons}_unsorted_cov.bdg

	fi
	awk -v OFS="\t" 'NR!=1{print $1,$2,$3,$4}' $out_dir/$name_cons/${name_cons}_unsorted_cov.bdg | sort -k1,1 -k2,2n > $out_dir/$name_cons/${name_cons}_cov.bdg
	rm $out_dir/$name_cons/${name_cons}_unsorted_cov.bdg
	$bdg2bwig $out_dir/$name_cons/${name_cons}_cov.bdg $sizeFile $out_dir/$name_cons/${name_cons}_cov.bw
	local i=0
	local files=()
	for dir in ${in_dir[@]}; # retrieving coverage values at maximum position to use as filter if necessary.
	do
		local name=${dir##*/}
		if [ $i -eq 0 ]; then
			bedtools intersect -a $out_dir/$name_cons/${name_cons}_narrow.bed -b $out_dir/$name/${name}_cov.bdg -wao | awk -v OFS="\t" '{print $1":"$2"|"$3,$7}' | sort -k1,1 -k2,2nr | sort -u -k1,1 > $out_dir/$name_cons/temp/tmp_max_${i}.txt
		else
			bedtools intersect -a $out_dir/$name_cons/${name_cons}_narrow.bed -b $out_dir/$name/${name}_cov.bdg -wao | awk -v OFS="\t" '{print $1":"$2"|"$3,$7}' | sort -k1,1 -k2,2nr | sort -u -k1,1 | awk '{print $2}' > $out_dir/$name_cons/temp/tmp_max_${i}.txt
		fi
		files+=("$out_dir/$name_cons/temp/tmp_max_${i}.txt")
		local i=$(($i+1))
	done
	paste ${files[@]} | sed "s/[|:]/\t/g" > $out_dir/$name_cons/${name_cons}_max.bed
	awk -v OFS="\t" '{mean=0;for(i=4;i<=NF;i++) {mean+=$i};print $1,$2,$3,mean/(i-3)}' $out_dir/$name_cons/${name_cons}_max.bed | sort -k4,4nr > $out_dir/$name_cons/${name_cons}_maxMean.bed
	rm -Rf $out_dir/$name_cons/*/
	LC_COLLATE=$tmpsave
}

#-------------------------------------------------------------------------------
peakcalling_MACS2 (){
	# peakcalling_MACS2 -id <PATH> -od <PATH> -g <INT> -s <INT> -phs <INT> -t <INT> -add "<STRING>" -kd <STRING> -bl <PATH> -size <PATH> -r <BOOL>
	# defining default parameters
	local Genome_length=120000000; local additionnal_args="";  local seedrandom=168159; local phs=200; local threads=1; local keepDup=1; local blacklist="/home/312.6-Flo_Re/312.6.1-Commun/data/A_thaliana_phytozome_v12/Greenscreen_19012023_merged.bed"; local sizeFile="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.chromsize"; local redo_analysis="false"
	# parsing arguments
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-id)
				local in_dir=$2
				echo "data directory set to: ${2}";shift 2;;
			-od)
				local out_dir=$2
				echo "output directory set to: ${2}";shift 2;;
			-g)
				local Genome_length=$2
				echo "Genome length (mappable) set to: ${2}";shift 2;;
			-s)
				local seedrandom=$2
				echo "random seed set to: ${2}"; shift 2;;
			-phs)
				local phs=$2
				echo "peaks half size set to: ${2}"; shift 2;;
			-t)
				local threads=$2
				echo "threads set to: ${2}"; shift 2;;
			-kd) 
				local keepDup=$2
				echo "level of duplication (auto, all or an integer as in MACS2 filterdup): ${2}"; shift 2;;
			-bl)
				local blacklist=$2
				echo "blacklist set to: ${2}"; shift 2;;
			-size)
				local sizeFile=$2
				echo "sizeFile set to: ${2}"; shift 2;;
			-r)
				local redo_analysis=$2
				echo "redo_analysis set to (true/false): ${2}"; shift 2;;
			-h)
				usage peakcalling_MACS2; return;;
			--help)
				usage peakcalling_MACS2; return;;
			*)
				echo "Error in arguments"
				echo $1; usage peakcalling_MACS2; return;;
		esac
	done
	# checking arguments for errors
	local Errors=0
	if [ -z $in_dir ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
	if [ -z $out_dir ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage peakcalling_MACS2; return 1; fi
	echo "Processing folder: ${in_dir##*/}" # treating replicate folder
	local out_dir=${out_dir}/${in_dir##*/}
	
	local folderdone=true
	# checking if replicate folder already has essential results files (Skipping the treatment of this folder if this is true
	if [ "$redo_analysis" == "true" ]; then 
		local folderdone=false
	elif [ ! -f ${out_dir}/${in_dir##*/}_max.bed ] || [ ! -f ${out_dir}/${in_dir##*/}_maxMean.bed ] || [ ! -f ${in_dir##*/}_cov.bdg.bed ]; then
		local folderdone=false
	fi
	mkdir -m 774 -p $out_dir
	# Saving Arguments
	echo $@ > $out_dir/Parameters_used.txt

	local controls=()
	local replicates=()
	# creating list of bam files for multiple controls and multiple replicates at once (usually only one replicate here as multiple replicate are dealt with mspc)
	echo "should print bam after"
	for bam in $(find $in_dir -name "*.filtered.sorted.bam")
	do
		echo "$bam"
		if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" == *"control"* ]]; then
			controls+=("$bam")
		else
			replicates+=("$bam")
		fi
	done
	local log=$out_dir/${in_dir##*/}.log;
	echo "" > $log
	local all_ok=true;
	
	if [ "$folderdone" == false ]; then
		mkdir -p -m 774 $out_dir/controls
		mkdir -p -m 774 $out_dir/replicates
		
		cd $out_dir/controls
		# checking if control files are present, not blocking the analysis if not, just changing macs2 launching parameters
		if [ ${#controls[@]} -eq 0 ]; then
			echo "No control files present"
			local all_ok=false;
		elif [ ${#controls[@]} -eq 1 ]; then
			echo "only one bam in control, no merge needed"
			cp ${controls[0]} control.bam
		else
			# pooling controls bam into a unique bam
			printf "%s\n" "${controls[@]}" > listbams.txt
			echo "concatenated controls: ${controls[@]}"
			$PATH_TO_SAMTOOLS/samtools merge --threads $(calc ${threads}-1) -b listbams.txt control.bam
		fi
		
		cd $out_dir/replicates
		# same check for replicates, but merging should not be usefull as we use MSPC to handle multiple replicates ( WARNING not tested with multiple replicates)
		if [ ${#replicates[@]} -eq 0 ]; then
			echo "No replicates present"
			local all_ok=false;
		elif [ ${#replicates[@]} -eq 1 ]; then
			echo "only one bam as replicate, no merge needed"
			cp ${replicates[0]} replicate.bam
		else
			echo "Replicates"
			printf "%s\n" "${replicates[@]}" > listbams.txt
			$PATH_TO_SAMTOOLS/samtools merge --threads $(calc ${threads}-1) -b listbams.txt replicate.bam
		fi
		R2_value="NA"
		# using minimal stats file of replicate to decide if analysis is paired ended or single ended as macs2 needs this info
		for bam in $(find $in_dir -name "*.filtered.sorted.bam")
		do
			if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
			R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${bam%.filtered.sorted.bam}.minimal.stats)
			fi
		done
		local formatBed="BED"
		if [ ${R2_value} != "NA" ]; then
			local formatBed="BEDPE" # format used by macs2
			# randsample is used to convert bam into bed format ( WARNING macs uses a special bed format, creating these files with macs2 suite is advised by macs2's creator)
			$PATH_TO_MACS/macs2 randsample -i $out_dir/replicates/replicate.bam -f BAMPE -p 100 -o $out_dir/replicates/replicate.bed
			if [ ${#controls[@]} -ne 0 ]; then
				$PATH_TO_MACS/macs2 randsample -i $out_dir/controls/control.bam -f BAMPE -p 100 -o $out_dir/controls/control.bed
			fi
		else
			$PATH_TO_MACS/macs2 randsample -i $out_dir/replicates/replicate.bam -f BAM -p 100 -o $out_dir/replicates/replicate.bed
			if [ ${#controls[@]} -ne 0 ]; then
				$PATH_TO_MACS/macs2 randsample -i $out_dir/controls/control.bam -f BAM -p 100 -o $out_dir/controls/control.bed
			fi
		fi
		cd $out_dir
		if [ $all_ok == true ]; then
			echo "Started MACS..."
			# call summit parameters is used by default as the summit of peaks are very usefull data for MSPC analysis afterwards
			echo "$PATH_TO_MACS/macs2 callpeak -t $out_dir/replicates/replicate.bed -c $out_dir/controls/control.bed -B -f $formatBed -n ${in_dir##*/} -g $Genome_length --call-summits --keep-dup $keepDup >> $log 2>&1;"
			$PATH_TO_MACS/macs2 callpeak -t $out_dir/replicates/replicate.bed -c $out_dir/controls/control.bed -B -f $formatBed -n ${in_dir##*/} -g $Genome_length --call-summits --keep-dup $keepDup >> $log 2>&1;
		elif [ ${#controls[@]} -eq 0 ]; then
			echo "Started MACS without controls"
			$PATH_TO_MACS/macs2 callpeak -t $out_dir/replicates/replicate.bed -B -f $formatBed -n ${in_dir##*/} -g $Genome_length --call-summits --keep-dup $keepDup --seed $seedrandom >> $log 2>&1;
		else
		   echo "Not all needed input present."
		   return 1
		fi
		# retrieving size of fragments used for normalisation of coverages
		if [ $formatBed == "BEDPE" ]; then
			fragment_length=$(cat $log | grep "fragment size = " | awk -v OFS="\t" '{print $12}')
		else
			fragment_length=$(cat $log | grep "predicted fragment length is" | awk -v OFS="\t" '{print $14}') # SE
		fi
		echo "=="
		echo "$fragment_length"
		echo "taking care of coverages values"
		# Using size of fragmet and bam file to recompute coverages along the genome.
		# We do not use the bdegraph generated by macs2 as they are normalize by size of smallest bam
		for bam in $(find $in_dir -name "*.filtered.sorted.dedup.bam")
		do
			if [[ ! -L $bam && -f $bam ]]; then
				if [ ${R2_value} == "NA" ]; then # checking for PE, or single ended reads
					if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
						#romain to dedup # I added ".dedup" in the for loop
						bedtools bamtobed -i $bam > $in_dir/${in_dir##*/}.bamtobed.bed
						awk -v fraglen=${fragment_length} -v OFS="\t" '$6=="+"{print $1,$2,$3+fraglen,$4,$5,$6}$6=="-"{print $1,$2-fraglen,$3,$4,$5,$6}' $in_dir/${in_dir##*/}.bamtobed.bed | awk -v OFS="\t" '$2<=0{print $1,1,$3,$4,$5,$6;next}{print $0}'  > $in_dir/${in_dir##*/}.ext.bed
						local libsize=$($PATH_TO_SAMTOOLS/samtools view -f 0 -c $bam)
						local scale=$(calc 1000000/$libsize)
						
						# sed -i 's/chr/Chr/g' $in_dir/${in_dir##*/}.ext.bed  ## Laura 27/05/2021
						tmpsave=$LC_COLLATE;LC_COLLATE=C
						bedtools genomecov -bga -scale $scale -i $in_dir/${in_dir##*/}.ext.bed -g $sizeFile | sed 's/Chr/chr/g' | sort -k1,1 -k2,2n > $out_dir/${in_dir##*/}_cov.bdg
						LC_COLLATE=$tmpsave
						$bdg2bwig $out_dir/${in_dir##*/}_cov.bdg $sizeFile $out_dir/${in_dir##*/}_cpm.bw
					fi
				else
					if [[ "$(echo "$bam" | tr '[:upper:]' '[:lower:]')" != *"control"* ]]; then
						local libsize=$($PATH_TO_SAMTOOLS/samtools view -f 0 -c $bam)
						local scale=$(calc 1000000/$libsize)
						#romain to dedup # I added ".dedup" in the for loop
						tmpsave=$LC_COLLATE;LC_COLLATE=C
						bedtools genomecov -bga -scale $scale -ibam $bam | sort -k1,1 -k2,2n > $out_dir/${in_dir##*/}_cov.bdg
						LC_COLLATE=$tmpsave
						$bdg2bwig $out_dir/${in_dir##*/}_cov.bdg $sizeFile $out_dir/${in_dir##*/}_cpm.bw
					fi
				fi
			fi
		done
		local nbpeaks=$(wc -l $out_dir/${in_dir##*/}_peaks.narrowPeak | awk '{print $1}')
		# checking if peaks are detected, as next computation doesn't like no peaks cases
		if [ $nbpeaks -eq 0 ]; then
			echo "ERROR no peaks in ${in_dir##*/}_peaks.narrowPeak created, MACS2 must have FAIL to compute peaks"
		elif [ ${#controls[@]} -eq 0 ]; then
			if [ $blacklist != "NA" ]; then
				echo "no controls but a blacklist has been submitted, filtering using these"
				bedtools intersect -v -a $out_dir/${in_dir##*/}_peaks.narrowPeak -b $blacklist -wa > $out_dir/${in_dir##*/}_filtered.narrowPeak
			else
				echo "no controls or blacklist submitted... filtering step omitted"
				cp $out_dir/${in_dir##*/}_peaks.narrowPeak $out_dir/${in_dir##*/}_filtered.narrowPeak
			fi
			local basename=${in_dir##*/}
		else
			echo "filtering out peaks detected in input"
			i=0
			for bdg in $(find $out_dir -name "*_treat_pileup.bdg")
			do
				# this whole section remove peaks that have 2 times mean genome coverage of input
				# this helps removing peaks that are present in the input and not only in the sample.
				i+=1
				local basename=${bdg%_treat_pileup.bdg}
				## Scaling factor for single-end data, counting every mapped read (bitwise flag = 0)
				ScaleTotalMappedReads=$(bc <<< "scale=6;1000000/$($PATH_TO_SAMTOOLS/samtools view -f 0 -c $out_dir/controls/control.bam)")
				tmpsave=$LC_COLLATE;LC_COLLATE=C
				bedtools genomecov -bga -ibam $out_dir/controls/control.bam -scale $ScaleTotalMappedReads | sort -k1,1 -k2,2n > $out_dir/${in_dir##*/}_control.bdg
				# Compute genome mean, weighted by length of each interval and doubled to secured the threshold
				local mean_cov=$(awk '{S+=$4*($3-$2);total+=$3-$2}END{printf "%.3f",(S/total)*10}' $out_dir/${in_dir##*/}_control.bdg)
				echo "tenfold mean coverage (threshold for forbidden regions) ${basename#${out_dir}/} $mean_cov"
				sort -k1,1 -k2,2n $out_dir/${in_dir##*/}_peaks.narrowPeak > $out_dir/${in_dir##*/}_sorted.narrowPeak
				bedtools intersect -a $out_dir/${in_dir##*/}_sorted.narrowPeak -b $out_dir/${in_dir##*/}_control.bdg -wb | awk -v OFS="\t" '{printf "%s\t%d\t%d\t%.2f\t%.3f\n",$11,$12,$13,$14,$9}' | sort -u -k1,1 -k2,2n | sort -k1,1 -k2,2n > $out_dir/${in_dir##*/}_controlonpeaks.bdg
				awk -v OFS="\t" -v Mean=${mean_cov} '$4 > Mean {print $0}' $out_dir/${in_dir##*/}_controlonpeaks.bdg > $out_dir/${in_dir##*/}_forbidden_regions.bed
				if [ $blacklist != "NA" ]; then 
					cat $blacklist >> $out_dir/${in_dir##*/}_forbidden_regions_unsorted.bed
					echo -e "\n" >> $out_dir/${in_dir##*/}_forbidden_regions_unsorted.bed
					awk -v OFS="\t" '{print $1,$2,$3}' $out_dir/${in_dir##*/}_forbidden_regions.bed >> $out_dir/${in_dir##*/}_forbidden_regions_unsorted.bed
					sort -k1,1 -k2,2n $out_dir/${in_dir##*/}_forbidden_regions_unsorted.bed > $out_dir/${in_dir##*/}_forbidden_regions.bed
				fi
				LC_COLLATE=$tmpsave
				$bdg2bwig $out_dir/${in_dir##*/}_control.bdg $sizeFile $out_dir/${in_dir##*/}_control.bw
				bedtools intersect -v -a $out_dir/${in_dir##*/}_peaks.narrowPeak -b $out_dir/${in_dir##*/}_forbidden_regions.bed -wa > $out_dir/${in_dir##*/}_filtered.narrowPeak
				## cleanup
				rm $out_dir/${in_dir##*/}_control.bdg
				rm $out_dir/${in_dir##*/}_controlonpeaks.bdg

			done
			if [ $i -eq 0 ]; then
				echo "ERROR no BEDGRAPH created, MACS2 must have FAIL to compute peaks"
			fi
		fi
		# removing temporary folders
		rm -R $out_dir/controls;
		rm -R $out_dir/replicates;
		# NOTE maybe we can check for multi-rep/single-rep here and make only one of these files per analysis ?
		awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' ${basename}_filtered.narrowPeak > ${basename}_peaks.bed
		awk -v OFS="\t" '{print $1, $2+$10-"'$phs'", $2+$10+"'$phs'"}' ${basename}_filtered.narrowPeak > ${basename}_narrow.bed
	else
		echo "Folder $out_dir already contains main results files. Skipping...";
	fi
}

# to be added in main_peakcalling if macs2 is used:
		# if [ $peakcaller == "MACS2" ]; then
		# 	echo "[INFO] - Using MACS2 for peakcalling"
		# 	# expected outputs of function peakcalling_MACS2 per replicate
		# 	echo "Launching MACS2 peakcalling" >> $out_dir/$name_cons/log.txt
		# 	list_peaks_mspc+=("$out_dir/$name/${name}_peaks.bed") 
		# 	list_bdg+=("$out_dir/$name/${name}_cov.bdg")
		# 	peakcalling_MACS2 -id $dir -od $out_dir -g $Genome_length -s $seedrandom -phs $phs -kd $keepDup -t $threads -size $sizeFile -r $redo_analysis -bl $blacklist
			
		# 	echo "Computing FRiP and other stats" >> $out_dir/$name_cons/log.txt
		# 	for elt in ${bam_ech[@]};
		# 	do
		# 		if [ ! -f $out_dir/$name/${name}_stats.txt ] || [[ $out_dir/$name/${name}_peaks.bed -nt $out_dir/$name/${name}_stats.txt ]]; then # checking if this step is already done to not do it again
		# 			if [[ $elt != *"control"* ]]; then # filter out all bam that contain control in his name.
		# 				# computing Fresquency of Reads in Peaks (FRiP)
		# 				total=$(samtools view -c $elt)
		# 				inpeak=$(bedtools sort -i $out_dir/$name/${name}_peaks.bed | bedtools merge -i stdin | bedtools intersect -u -a $elt -b stdin -ubam | samtools view -c)
		# 				FRIP=$(calc $inpeak/$total*100)
		# 				echo "$name $FRIP% $inpeak / $total" | tee -a $out_dir/$name_cons/log.txt
		# 				echo "$name $FRIP% $inpeak / $total" > $out_dir/$name/${name}_FreqReadInPeak.txt
		# 				R2_value=$(awk -v FS=" " 'NR>1{print $2}' ${elt%.filtered.sorted.bam}.minimal.stats)
		# 				if [ $R2_value == "NA" ]; then # number of tags/fragsment, way of rzetrienving this info is PE/SE dependant
		# 					totTags=$(grep "total tags in treatment" $out_dir/$name/${name}.log | cut -d " " -f 15)
		# 					if [ $keepDup == "all" ]; then
		# 						filtTags=$totTags
		# 					else
		# 						filtTags=$(grep "tags after filtering in treatment" $out_dir/$name/${name}.log | cut -d " " -f 16)
		# 					fi
		# 				fi
		# 				if [ $R2_value != "NA" ]; then
		# 					totTags=$(grep "total fragments in treatment" $out_dir/$name/${name}.log | cut -d " " -f 15)
		# 					if [ $keepDup == "all" ]; then
		# 						filtTags=$totTags
		# 					else
		# 						filtTags=$(grep "fragments after filtering in treatment" $out_dir/$name/${name}.log | cut -d " " -f 16)
		# 					fi
		# 				fi
		# 				filtPeaks=$(wc -l $out_dir/$name/${name}_filtered.narrowPeak | cut -d " " -f 1)
		# 				echo "Sample totalTags filtTags filtPeaks FRIP" > $out_dir/$name/${name}_stats.txt
		# 				echo $name $totTags $filtTags $filtPeaks $FRIP >> $out_dir/$name/${name}_stats.txt
		# 			fi
		# 		fi
		# 	done
		# fi




#-------------------------------------------------------------------------------
comparison(){
	# comparison -n -od -id -f -he
	while  [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-n)
				local names=("${!2}")
				echo "name of directory for dataset 1 is: ${names[@]}";shift 2;;
			-od)
				local result=$2
				echo "output directory set to: ${2}";shift 2;;
			-id)
				local data=$2
				echo "general data directory set to: ${2}";shift 2;;
			-f)
				local filterCov=$2
				echo "coverage must be at least >${2} in all samples to be considered as peaks";shift 2;;
			-he)
				local filterHeight=$2
				echo "height must be at least >${2} in all samples to be considered as peaks";shift 2;;
			-h)
				usage comparison; return;;
			--help)
				usage comparison; return;;
			*)
				echo "Error in arguments"
				echo $1; usage comparison; return;;
		esac
	done
	local Errors=0
	if [ -z $names ]; then echo "ERROR: -n argument needed or need to be a LIST"; Errors+=1; fi
	if [ -z $result ]; then echo "ERROR: -od argument needed"; Errors+=1; fi
	if [ -z $data ]; then echo "ERROR: -id argument needed"; Errors+=1; fi
	if [ -z $filterCov ]; then echo "no filter on coverage applied"; local filterCov=0; fi
	if [ -z $filterHeight ]; then echo "no filter on peak height applied"; local filterHeight=0; fi
	if [ $Errors -gt 0 ]; then usage comparison; return 1; fi


	local out_dir=$result ## LT 09/01; was $result/ before
	
	mkdir -p -m 774 ${out_dir}
	local peaks=()
	echo "here are the names: ${names[@]}" ## LT 09/01
	for ((i=0;i<${#names[@]};i++)) # get all bdgs & peaks files in two lists
	do
		covs+=("$data/${names[i]}/${names[i]}_cov.bdg")
		peaks+=("${out_dir}/${names[i]}_peaks.bed")
		if [ $filterHeight -eq 0 ];then
			echo ${names[i]} ## LT 09/01
			awk -v OFS="\t" -v name=${names[i]} '{print $1,$2,$3,name}' $data/${names[i]}/${names[i]}_narrow.bed  > ${out_dir}/${names[i]}_peaks.bed
		else # we filter by height of maximum of peaks
			awk -v filt=$filterHeight -v OFS="\t" -v keep=1 '{for(i=4;i<=NF;i++) {if($i<=filt){keep=0}};if(keep==1){print $0};keep=1}' 	$data/${names[i]}/${names[i]}_max.bed > $data/${names[i]}/${names[i]}_narrow.heightFiltered.bed
			awk -v OFS="\t" -v name=${names[i]} '{print $1,$2,$3,name}' $data/${names[i]}/${names[i]}_narrow.heightFiltered.bed  > ${out_dir}/${names[i]}_peaks.bed
		fi
	done

	local str_peaks=`join_by " " "${peaks[@]}"`


	awk '{print}' "${peaks[@]}" | sort -k1,1 -k2,2n > ${out_dir}/all_peaks.bed ## LT 09/01
	# cat ${str_peaks} | sort -k1,1 -k2,2n > ${out_dir}/all_peaks.bed ## commented LT 09/01

	python3 ${merge_peaks_Nsets} -b ${out_dir}/all_peaks.bed -o ${out_dir}/merged.bed
	echo "peaks merged"

	local inter_cov=("${out_dir}/merged.bed")
	local header=("chr" "begin" "end" "name")
	for ((i=0;i<${#names[@]};i++))
	do
		bedtools intersect -a ${out_dir}/merged.bed -b ${covs[i]} -wa -wb -sorted -loj | awk  -v OFS="\t" '$6 != "-1" {print $0} $6=="-1" {print $1,$2,$3,$4,$1,$2,$3,0}' > ${out_dir}/${names[i]}.inter
		python3 $compute_coverage -i ${out_dir}/${names[i]}.inter -m
		Noise=$(awk '{sum+=$4*($3-$2);total+=$3-$2}END{print sum/total}' ${covs[i]})
	# 	awk -v noise=$Noise '{print ($5/($3-$2))-noise}' ${out_dir}/${names[i]}.inter.cov | awk '$1>=0.0{print $0;next}{print 0.0}' > ${out_dir}/${names[i]}.inter.cov.normed
		awk -v noise=$Noise '{print ((1000*$5)/($3-$2))}' ${out_dir}/${names[i]}.inter.cov | awk '$1>=0.0{print $0;next}{print 0.0}' > ${out_dir}/${names[i]}.inter.cov.normed
		inter_cov+=("${out_dir}/${names[i]}.inter.cov.normed")
		header+=("${names[i]}")
		echo "${names[i]} treated"
	done
	# echo "${inter_cov[@]}" ## LT 09/01
	local str_inter_cov=`join_by " " "${inter_cov[@]}"`
	local str_header=`join_by "	" "${header[@]}"`
	echo ${str_header}
	if [ $filterCov -eq 0 ];then
		paste ${inter_cov[@]} | sed "1i${str_header}" >  ${out_dir}/table_peaks.csv ## LT 10/01/2022
	# 	paste ${str_inter_cov} | sed "1i${str_header}" >  ${out_dir}/table_peaks.csv
	else
		paste ${str_inter_cov} | awk -v filter=$filterCov -v OFS="\t" '$5>=filter || $6>=filter {print $0}' |  sed "1ichr\tbegin\tend\tname\t$name1\t$name2" | sed "1i${str_header}" >  ${out_dir}/table_peaks.csv
	fi
	awk -v OFS="\t" 'NR==1{print $0,"ratio"$5,"ratio"$6,"ratio"$7;next}{print $0,$5/$5,$6/$5,$7/$5}' ${out_dir}/table_peaks.csv > ${out_dir}/table_peaks_ratios.csv

}
# comparison) 
# echo -e "
# ==========
# DEPRECATED SINCE 2023
# usage: comparison -n <LIST of STRING> -od <PATH> -id <PATH> 
#        -f [INT] -he [FLOAT] 

# general infos: comparison of peaks coverage from any N samples or replicats,
#                filtered (optional) by the height of maximum or global peak 
#                coverages.

# -- Mandatory arguments:
#     -n      STRING  :    list of names of directory for dataset 1
#     -od     PATH    :    Set the Output directory, where you will find the 
#                          BEDGRAPH and NarrowPeak files created by this 
#                          programs.
#     -id     PATH    :    Set the directory where the BED/BEDGRAPH files are
#                          stored.

# -- optional arguments :
#     -f      INT     :    Threshold of coverage -in RPKM-. Filters out peaks
#                          for which for both samples do not pass the threshold.
#     -he     INT     :    Threshold on the maximum height of coverage -in RPKM-
#                          Filters out peaks for which both samples have maximum
#                          of coverage below this value.
# \n";;



# Deprecated by JL the 25/03/2026 for a more performant version.
compute_distribution(){
	# FUNCTION: compute_distribution
	# DESCRIPTION:
	#   This function calculates and visualizes the distribution of motif scores (TFFM or PWM)
	#   across a genome. It processes each chromosome in parallel to generate scores,
	#   samples a subset of these scores, and creates an interactive distribution plot.
	#   This is particularly useful for analyzing the global occurrence of motifs in a genome 
	#   and determining threshold for downstream analysis.
	#
	# USAGE:
	#   compute_distribution -m <file> -n <STRING> -od <PATH> -g [FILE] -r [INT]
	#         -c [LIST of STRING]
	#
	# ARGUMENTS:
	#   -m          : Path to the matrix file (TFFM in `.xml` format or PWM in `.pfm` format). Required.
	#   -od         : Directory where results will be saved. Required.
	#   -n          : Prefix for output files and directories. Required.
	#   -g          : Path to the genome FASTA file. Optional. Default: "/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas".
	#   -r          : Identifier for the interactive plot request. Optional. Default: 0.
	#   -c          : List of chromosomes to analyze. Optional. Default: ("chr1" "chr2" "chr3" "chr4" "chr5").
	#   -h, --help  : Display this help message.
	#
	# DEPENDENCIES 
		# - ../bin/scores.py 
		# - get_all_score_tffm.py
		# - Environment with pyhon version 2.7 for TFFM
		# - ../bin/InteractiveDistrib.r
	# NOTES:
	#   - For TFFM matrices, the function computes scores for the entire genome at once.
	#   - For PWM matrices, the function processes each chromosome in parallel to speed up computation.
	#   - The function samples 10% of the scores from each chromosome to create a manageable dataset.

	echo "
	____ ____ _  _ ___  _  _ ___ ____     ___  _ ____ ___ ____ _ ___ 
	|    |  | |\/| |__] |  |  |  |___     |  \ | [__   |  |__/ | |__]
	|___ |__| |  | |    |__|  |  |___ ___ |__/ | ___]  |  |  \ | |__]

	"
	local genome="/home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; 
	local chrset=("chr1" "chr2" "chr3" "chr4" "chr5"); 
	local request=0
	local Argline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-g)
				local genome=$2
				echo "-> FASTA of the genome set to: ${2}";shift 2;;
			-m)
				local matrice=$2
				echo "-> Matrix files set to: ${matrice}";shift 2;;
			-od)
				local results=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-n)
				local name=$2
				echo "-> Names associated set to: ${name}";shift 2;;
			-r)
				local request=$2;shift 2;;
			-c)
				local chrset=("${!2}"); echo "-> Chromosomes in genome: ${!2}"
				shift 2;;
			-h)
				echo "you asked for help" ;Errors+=1; shift 1;;
			--help)
				echo "you asked for help" ;Errors+=1; shift 1;;
			*)
				echo "Error in arguments"
				echo $1; usage compute_distribution ; break;;
		esac
	done

	# checking arguments
	local Errors=0
	if [ -z $genome ]; then echo "-g argument not used, assuming A.thaliana is used: /home/312.6-Flo_Re/312.6.1-Commun/data/tair10.fas"; fi
	if [ $Errors -gt 0 ]; then usage compute_space; return 1; fi
	
	mkdir -p -m 774 $results/${name} $results/${name}/scores/chr1 $results/${name}/scores/chr2 $results/${name}/scores/chr3 $results/${name}/scores/chr4 $results/${name}/scores/chr5

	# scores_prog=/home/312.6-Flo_Re/312.6.1-Commun/Jeremy/Opti_scores.py
	# Computing scores
	if [[ $matrice == *".xml"* ]]; then # TFFM 
		$Python_TFFM $tffm_all_scores -o ${results}/${name}/scores/${name}_tffm_scores_pos.tsv -pos $genome -t ${matrice}
		local scores=${results}/${name}/scores/${name}_tffm_scores_pos.tsv
	elif [[ $matrice == *".pfm"* ]]; then # PWM
		local scores=${results}/${name}/scores/$(basename $genome).scores
		if [[ ! -f $scores ]] || [[ ${matrice} -nt $scores ]] ; then
			python3 ${scores_prog} -m ${matrice} -f $genome -o ${results}/${name}/scores/chr1 -chr ${chrset[0]} &
			python3 ${scores_prog} -m ${matrice} -f $genome -o ${results}/${name}/scores/chr2 -chr ${chrset[1]} &
			python3 ${scores_prog} -m ${matrice} -f $genome -o ${results}/${name}/scores/chr3 -chr ${chrset[2]} &
			python3 ${scores_prog} -m ${matrice} -f $genome -o ${results}/${name}/scores/chr4 -chr ${chrset[3]} &
			python3 ${scores_prog} -m ${matrice} -f $genome -o ${results}/${name}/scores/chr5 -chr ${chrset[4]} &
			wait
			
			# merge all files together
			awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print $0}' ${results}/${name}/scores/chr1/$(basename $genome).scores > ${results}/${name}/scores/$(basename $genome).scores
			awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print $0}' ${results}/${name}/scores/chr2/$(basename $genome).scores >> ${results}/${name}/scores/$(basename $genome).scores
			awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print $0}' ${results}/${name}/scores/chr3/$(basename $genome).scores >> ${results}/${name}/scores/$(basename $genome).scores
			awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print $0}' ${results}/${name}/scores/chr4/$(basename $genome).scores >> ${results}/${name}/scores/$(basename $genome).scores
			awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .1 || FNR==1) print $0}' ${results}/${name}/scores/chr5/$(basename $genome).scores >> ${results}/${name}/scores/$(basename $genome).scores

			rm -Rf $results/${name}/scores/chr1 $results/${name}/scores/chr2 $results/${name}/scores/chr3 $results/${name}/scores/chr4 $results/${name}/scores/chr5
		fi
	fi
	/home/prog/R/R-4.4.0/el8/bin/Rscript $interactivedistrib $scores $results/${name}/ $request
}


#-------------------------------------------------------------------------------
replicates_comparisons (){
	# FUNCTION: replicates_comparisons
	# DESCRIPTION:
	#   This function compares replicates for peak calling and coverage computation. 
	#   It supports two normalization methods: "Not Input Normalized" and "Input Re-Normalized". 
	#
	# USAGE:
	#   replicates_comparisons -p <PATH> -nb <INT> -od <PATH>
	#                          -gd <PATH> -n <STRING> -c [STRING]
	#
	# ARGUMENTS:
	#   -p          : Path to the consensus peaks file. Required.
	#   -nb         : Number of replicates in the sample. Required.
	#   -od         : Output directory for results. Required.
	#   -gd         : General directory containing `Peakcalling` and `Mapping` subdirectories. Required.
	#   -n          : Name of the sample. Required.
	#   -c          : Color for dots in plots (hex code, e.g., `#FF0000` or `FF0000`). Optional. Default: `#000000`.
	#   -h, --help  : Display usage information for this function.
	# 
	# DEPENDENCIES:
	#   - compute_rpkmrip_rpkmril
	#   - Pairwise_comps (R script)

	echo "
	____ ____ ___      ____ ____ _  _ ___  ____ ____ _ ____ ____ _  _
	|__/ |___ |__]     |    |  | |\/| |__] |__| |__/ | [__  |  | |\ |
	|  \ |___ |    ___ |___ |__| |  | |    |  | |  \ | ___] |__| | \|
	"
	#TODO: add option to show consensus peaks and/or all peaks from all replicates. (JL 09/2025; from comments of FP) ## JL: I'll do a branch for this soon
	local color="#000000";
	local Argsline=$@
	while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
		case $1 in
			-p)
				local peaks=$2
				echo "-> Consensus peaks file set to: ${2}";shift 2;;
			-nb)
				local nbrep=$2
				echo "-> Number of replicates in the sample set to: ${2}";shift 2;;
			-od)
				local result=$2
				echo "-> Output directory set to: ${2}";shift 2;;
			-gd)
				local general_dir=$2
				echo "-> General result directory set to: ${2}";shift 2;;
			-n)
				local name=$2
				echo "-> Name of sample set to: ${2}";shift 2;;
			-c)
				local color=$2;
				echo "-> Color set for dots: $color"; shift 2;;
			-h)
				usage replicates_comparisons ; return;;
			--help)
				usage replicates_comparisons ; return;;
			*)
				echo "[ERROR] - In argument: $1 $2"
				usage replicates_comparisons ; return;;
		esac
	done
	local Errors=0
	if [ -z $peaks ]; then echo "ERROR: -p argument is needed"; Errors+=1; fi
	if [ -z $nbrep ]; then echo "ERROR: -nb argument is needed"; Errors+=1; fi
	if [ -z $result ]; then echo "ERROR: -od argument is needed"; Errors+=1; fi
	if [ -z $name ]; then echo "ERROR: -n argument is needed"; Errors+=1; fi
	if [ -z $general_dir ]; then echo "ERROR: -gd argument is needed"; Errors+=1; fi
	if [ $Errors -gt 0 ]; then usage replicates_comparisons ; return 1; fi

	[[ ${color} =~ ^#.* ]] || color="#${color}"

	# initialize lists
	local rep_list=()
	local bdg_list=()
	local remove_list=()
	mkdir -p -m 774 $result/$name/NotInputNormalized $result/$name/InputReNormalized
	local log=$result/$name/log.txt

	echo "" > $log
	echo $Argsline >> $log

	# prepare lists of replicates files
	for rep in $(seq $nbrep)
	do
		rep_list+=("${name}rep${rep}")
		bdg_list+=("$general_dir/Peakcalling/${name}rep${rep}/${name}rep${rep}_cpm.bdg")
		rep_peaks_list+=("$general_dir/Peakcalling/${name}rep${rep}/${name}rep${rep}_narrow.bed")
		remove_list+=("$result/$name/NotInputNormalized/${name}rep${rep}_filt.cov.bed" "$result/$name/InputReNormalized/${name}rep${rep}_filt.cov.bed" "$result/$name/NotInputNormalized/${name}rep${rep}_cpmrip.bdg" "$result/$name/InputReNormalized/${name}rep${rep}_cpmrip.bdg")
	done

	#concatenate peaks from all replicates, keeping track of their origin
	list_data=("$general_dir/Peakcalling")
	list_bamdir=("$general_dir/Mapping")

	# Not Input Normalized method
	echo "[INFO] - Computing coverages with RPKM-RiP method (not Input normalized)" | tee -a $log
	compute_rpkmrip_rpkmril -p $peaks -bd list_bamdir[@] -pd list_data[@] -sn rep_list[@] -o $result/$name/NotInputNormalized/ -m "inPeaks"
	$PATHRscript $Pairwise_comps $result/$name/NotInputNormalized/ $result/$name/NotInputNormalized/peaks_perSample_rpkminPeaks.txt $color

	# Input ReNormalized method
	echo "[INFO] - Computing coverages with RPKM-RiL method (Input Re-Normalized)" | tee -a $log
	compute_rpkmrip_rpkmril -p $peaks -bd list_bamdir[@] -pd list_data[@] -sn rep_list[@] -o $result/$name/InputReNormalized/ -m "inLibs"
	$PATHRscript $Pairwise_comps $result/$name/InputReNormalized/ $result/$name/InputReNormalized/peaks_perSample_rpkminLibs.txt  $color

	# Cleaning 
	echo "[INFO] - Cleaning temporary files" | tee -a $log
	for f in ${remove_list[@]}; do
		if [ -f $f ]; then
			rm $f
		else
			echo "No $f to remove !" | tee -a $log
		fi
	done
	# rm ${remove_list[@]}
	rm $result/$name/*/tmp_peaks
	rm $result/$name/NotInputNormalized/tmpFiltTags.txt || echo "No $result/$name/NotInputNormalized/tmpFiltTags.txt to remove !"
	rm $result/$name/InputReNormalized/tmpFiltTags.txt || echo "No $result/$name/InputReNormalized/tmpFiltTags.txt to remove !"
}
