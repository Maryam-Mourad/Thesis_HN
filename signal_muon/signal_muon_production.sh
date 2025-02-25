IDs=(41 42 43 44 434)
masses=(2.2 3.0 4.0 5.0 4.5)
#ctaus=(3.041 0.645 0.153 0.050)
#ctaus=(4628.421 879.474 161.646 41.739) 
ctaus=(1131 219 46.4 14.18 24.8) #coupling square times ctau times 10E5 [mm] #baseline_ctau

nEvents=1000000
generation=0
if [ $1 ]; then
    nEvents=$1
fi

if [ $2 ]; then
    generation=$2
fi





let temp=${#masses[@]}-1
for massPoint in $(seq 0 $temp)
do
    echo ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    echo Mass point to use: ${masses[$massPoint]}
    echo ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    echo 
    echo 
    if [ -r HN${masses[$massPoint]} ]; then
        echo HN directory exists.
        cp mastercode_signal_muon.py HN${masses[$massPoint]}/
    else
        mkdir HN${masses[$massPoint]}
        echo HN directory created.
        mkdir HN${masses[$massPoint]}/plots
        cp signal_muon.config HN${masses[$massPoint]}/signal_muon_mass${masses[$massPoint]}.config
        cp signal_muon.decay HN${masses[$massPoint]}/signal_muon_mass${masses[$massPoint]}.decay
        cp mastercode_signal_muon.py HN${masses[$massPoint]}/
        python3 editDecay.py ${masses[$massPoint]}
    fi
    python3 addHN.py ${IDs[$massPoint]} ${masses[$massPoint]} ${ctaus[$massPoint]}
    cd HN${masses[$massPoint]}
    echo 
    echo 
    if [ $generation = generate ]; then
        echo generating events ...
        $RAPIDSIM_ROOT/build/src/RapidSim.exe $RAPIDSIM_ROOT/HN_simulation/signal_muon/HN${masses[$massPoint]}/signal_muon_mass${masses[$massPoint]} $nEvents 1
        mv signal_muon_mass${masses[$massPoint]}_tree.root /eos/user/m/mmourad/HN_simulation/trees/
    else
        echo Skipping MC generation and moving to plotting ...
    fi
    echo
    echo
    echo ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    echo running the mastercode for this mass point ...
    #python3 mastercode_signal_muon.py ${masses[$massPoint]} $massPoint $nEvents ${ctaus[$massPoint]}
    echo finished mass value ${masses[$massPoint]}
    echo ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    echo 
    echo 
    cd ..
done

hadd -f histograms.root HN*/histograms_mass_*.root

if [ -r plots ]; then
    echo
else
    mkdir plots
    echo plots directory created.
fi

python3 plotter.py ${masses[@]}


hadd -f histograms_muon.root histograms.root ../gamma2/histograms.root ../gamma42/histograms.root ../gamma9/histograms.root
cd ..

if [ -r plots_muon ]; then
    echo
else
    mkdir plots_muon
    echo plots directory created.
fi


python3 plotter_muon.py ${masses[@]}


echo All done ":)"