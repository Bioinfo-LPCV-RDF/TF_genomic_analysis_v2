bash | find -name AUC.txt | while read -r line; do grep "PWM" $line; done
