



import argparse #v1.1
from pptx import Presentation #v0.6.21
from pptx.util import Inches, Pt, Cm



parser = argparse.ArgumentParser()

parser.add_argument("--motifPWM", "-pwm", help='')
parser.add_argument("--numberpeaks", "-pknb", help='')
parser.add_argument("--ROCs", "-roc", help='')
parser.add_argument("--replicateComparison", "-rep", help='', default="NA")
parser.add_argument("--spacing", "-spa", help='')
parser.add_argument("--samplename", "-n", help='')
parser.add_argument("--out", "-o", help='path and filename without extension')


args = parser.parse_args()

def main(motifPWM, ROCs, replicateComparison, spacing, samplename, out, numberpeaks):
	# Create a new PowerPoint presentation
	prs = Presentation()

	# Add a slide to the presentation
	slide = prs.slides.add_slide(prs.slide_layouts[6])
	
	heightREP = Cm(8.44)
	topPWM = Cm(8); heightPWM = Cm(4)
	topROC = Cm(12); heightROC = Cm(7)
	leftSPA = Cm(12); topSPA = Cm(7); heightSPA = Cm(12)
	
	if replicateComparison != "NA":
		picREP = slide.shapes.add_picture(replicateComparison, left=0, top=0, height=heightREP)
	else:
		textbox = slide.shapes.add_textbox(Cm(2), Cm(4), 200, 100)
		textbox.text = "No replicates for this sample"
	picPWM = slide.shapes.add_picture(motifPWM, left=0, top=topPWM, height=heightPWM)
	picROC = slide.shapes.add_picture(ROCs, left=0, top=topROC,height=heightROC)
	picSPA = slide.shapes.add_picture(spacing, left=leftSPA, top=topSPA, height=heightSPA)

	## Add a textbox to the slide
	textbox = slide.shapes.add_textbox(Cm(10), 0, 200, 100)

	## Add text to the textbox
	textbox.text = "Basic Analysis: "+samplename+"\n Number of Consensus peaks: "+str(numberpeaks)
	

	# Save the PowerPoint presentation
	prs.save(out)




main(args.motifPWM, args.ROCs, args.replicateComparison, args.spacing, args.samplename, args.out, args.numberpeaks)