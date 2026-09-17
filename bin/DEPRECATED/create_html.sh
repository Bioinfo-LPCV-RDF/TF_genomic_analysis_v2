#!/bin/bash




create_html(){
local getCov="yes"; 
while [ $# -ge 1 ] && [[ -n $1 ]] && [[ $1 != "\n" ]] ; do
	case $1 in
		-n)
			local NameTF=$2
			echo "TFname set to ${2}";shift 2;;
		-od)
			local out_dir=$2
			echo "output directory set to: ${2}"; shift 2;;
		--motif-directory)
			local motif_dir=$2
			echo "directory where logo can be found set to: ${2}"; shift 2;;
		--reps-comp-directory)
			local reps_dir=$2
			echo "directory for comparisons of repliactes se to: ${2}"; shift 2;;
		--rocs-directory)
			local roc_dir=$2
			echo "directory where ROCs can be found set to: ${2}"; shift 2;;
		-h)
			usage create_html; return;;
		--help)
			usage create_html; return;;
		*)
			echo "Error in arguments"
			echo $1; usage create_html; return;;
	esac	
done
local Errors=0
if [ -z $NameTF ]; then echo "ERROR: -n argument needed"; Errors+=1; fi

if [ $Errors -gt 0 ]; then usage create_html; return 1; fi	

mkdir -p -m 774 $out_dir
mkdir -p -m 774 $out_dir/$NameTF


#### IDEES
# mieux vaut copier les images dans un dossier attenant au html résultat pour éviter que le html soit inutilisable si l'individu ne possède pas les droit de lecture sur les dossiers indiqués



#### BEGIN Treatment for motif
logo_list=()
logoRC_list=()
info_list=()
continue_search=true;found=false; rcfound=false
i=1
while $continue_search; do
	if [ -f $motif_dir/meme/meme_out/logo${i}.png ]; then 
# 		logo_list+=("$motif_dir/meme/meme_out/logo${i}.png")
		cp $motif_dir/meme/meme_out/logo${i}.png $out_dir/$NameTF/logo${i}.png
		logo_list+=("$NameTF/logo${i}.png")
		found=true
	fi
	if [ -f $motif_dir/meme/meme_out/logo_rc${i}.png ]; then 
# 		logoRC_list+=("$motif_dir/meme/meme_out/logo_rc${i}.png")
		cp $motif_dir/meme/meme_out/logo_rc${i}.png $out_dir/$NameTF/logoRC${i}.png
		logoRC_list+=("$NameTF/logoRC${i}.png")
		rcfound=true 
	fi
	if $found && $rcfound; then
		totalseq=$(wc -l $motif_dir/meme/seqs-centered | awk '{print $1/2}')
		info_list+=("$(grep 'MOTIF' $motif_dir/meme/meme_out/meme.txt | grep -v 'SEQUENCE\|BL\|SUMMARY' | grep "MEME-${i}" | awk -v OFS="</p><p>" -v total=$totalseq '{print "<p>number of detected sites: "$9"/"total,"Log Likelihood Ratio: "$12,"E-value: "$15"</p>"}')")
	else
		continue_search=false
		
	fi
	i=$(($i+1))
	found=false
	rcfoun=false
done

table_motif="<table><tbody>"
for (( c=0; c<${#info_list[@]}; c++ ))
do 
	if [ $c -eq 0 ]; then
		table_motif+="<tr class=\"selected\">"
	else
		table_motif+="<tr>"
	fi
	table_motif+="<td><img src=\"${logo_list[$c]}\" width=\"200\" height=\"120\"/></td>"
	table_motif+="<td><img src=\"${logoRC_list[$c]}\" width=\"200\" height=\"120\"/></td>"
	table_motif+="<td><p>Motif MEME-$(($c+1))</p>${info_list[$c]}</td>"
	table_motif+="</tr>"
done
table_motif+="</tbody></table>"

#### END Treatment for motif

#### BEGIN Treatment for replicate comparison
cp $reps_dir/IRN/between_samples_logScale.png $out_dir/$NameTF/IRN_repcomp.png
cp $reps_dir/NIN/between_samples_logScale.png $out_dir/$NameTF/NIN_repcomp.png
#### END Treatment for replicate comparison

#### BEGIN Treatment for ROCs
cp $roc_dir/ROC.png $out_dir/$NameTF/ROC.png
#### END

#### BEGIN Treatment for 

#### END

#### BEGIN Treatment for 

#### END



touch $out_dir/$NameTF.html

cat > $out_dir/$NameTF.html << EOF

<!DOCTYPE html>
<html>
	<head>
		<link rel="stylesheet" href="$NameTF.css">
		<title>${NameTF} TF analysis</title>
		<script src="https://code.jquery.com/jquery-3.6.0.min.js" integrity="sha256-/xUj+3OJU5yExlq6GSYGSHk7tPXikynS7ogEvDej/m4=" crossorigin="anonymous"></script>
		<script> 
			function changeReps()
			{
			var img = document.getElementById("Reps_comp_image");
			img.src="$NameTF/NIN_repcomp.png";
			document.getElementById("Reps_comp_txt").innerHTML="RPKM-RiP,<br> Reads per kilo base per million mapped reads<br>Only Reads in Peaks are considered";
			return false;
			}
			function changeRepsBack()
			{
			var img = document.getElementById("Reps_comp_image");
			img.src="$NameTF/IRN_repcomp.png";
			document.getElementById("Reps_comp_txt").innerHTML="RPKM-RiL,<br> Reads per kilo base per million mapped reads<br>All Reads in Library are considered";
			return false;
			}
		</script>
	</head>
	<body>
		<h1 align="center">${NameTF}</h1>
		<table align="center">
			<tbody>
				<tr>
					<td>
						$table_motif
					</td>
					<td>
						<table>
							<tbody>
								<tr>
									<td>
									<p align="left" id="Reps_comp_txt">RPKM-RiL,<br> Reads per kilo base per million mapped reads<br>All Reads in Library are considered</p>
									</td>
									<td align="right">
										<button id="clickme" onclick="changeReps();">RPKM-RiP</button>
										<button id="clickme" onclick="changeRepsBack();">RPKM-RiL</button>
									</td>
								</tr>
								<tr>
									<td><img id="Reps_comp_image" src="$NameTF/IRN_repcomp.png" width=500 height=500/>
									</td>
								</tr>
							</tbody>
						</table>
					</td>
				<tr>
			</tbody>
		</table>
		<table>
			<tbody>
				<tr>
					<td><img src="$NameTF/ROC.png" width=500 height=500>
					<td>
				</tr>
			</tbody>
		</table>
	</body>
</html>

EOF
#     <iframe width="660" height="425" src="${htmlimagevar}" seamless></iframe>
#     <p>This is from index.html</p>
#     <img class="fit-picture"
#      src="${imagevar}"
#      alt="Grapefruit slice atop a pile of other slices">
#      <p>this is the explanation MF</p>



# <script> 
#     \$(function(){
#       \$("#includedContent").load("${htmlimagevar}"); 
#     });
#     </script> 

touch $out_dir/$NameTF.css
cat > $out_dir/$NameTF.css << EOF
.selected{
background-color: rgba(0, 0, 0, .1);
}

EOF
}



create_html -n ARF2FL \
			-od /home/312.3-StrucDev/312.3.1-Commun/ARF-anr/DAP_052022/ARF/ARF5/Summary \
			--motif-directory /home/312.3-StrucDev/312.3.1-Commun/ARF-anr/DAP_052022/ARF/ARF5/Motifs/ARF2FL \
			--reps-comp-directory /home/312.3-StrucDev/312.3.1-Commun/ARF-anr/DAP_052022/ARF/ARF5/Reps_comparison/ARF2FL \
			--rocs-directory /home/312.3-StrucDev/312.3.1-Commun/ARF-anr/DAP_052022/ARF/ARF5/ROCs/ARF2FL