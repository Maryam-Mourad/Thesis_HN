IDs=(41 42 43 44)
masses=(2.2 3.0 4.0 5.0)
ctaus=(3.041 0.645 0.153 0.050)

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
        cp mastercode_signal_tau.py HN${masses[$massPoint]}/
    else
        mkdir HN${masses[$massPoint]}
        echo HN directory created.
        mkdir HN${masses[$massPoint]}/plots
        cp signal_tau.config HN${masses[$massPoint]}/signal_tau_mass${masses[$massPoint]}.config
        cp signal_tau.decay HN${masses[$massPoint]}/signal_tau_mass${masses[$massPoint]}.decay
        cp mastercode_signal_tau.py HN${masses[$massPoint]}/
        python3 editDecay.py ${masses[$massPoint]}
    fi
    python3 addHN.py ${IDs[$massPoint]} ${masses[$massPoint]} ${ctaus[$massPoint]}
    cd HN${masses[$massPoint]}
    echo 
    echo 
    if [ $generation = generate ]; then
        echo generating events ...
        $RAPIDSIM_ROOT/build/src/RapidSim.exe $RAPIDSIM_ROOT/HN_simulation/signal_tau/HN${masses[$massPoint]}/signal_tau_mass${masses[$massPoint]} $nEvents 1
        mv signal_tau_mass${masses[$massPoint]}_tree.root /eos/user/m/mmourad/HN_simulation/trees/
    else
        echo Skipping MC generation and moving to plotting ...
    fi
    echo
    echo
    echo ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    echo running the mastercode for this mass point ...
    python3 mastercode_signal_tau.py ${masses[$massPoint]} $massPoint $nEvents
    echo finished mass value ${masses[$massPoint]}
    echo ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    echo 
    echo 
    cd ..
done

hadd -f histograms.root HN*/histograms_mass_*.root
python3 plotter.py ${masses[@]}

hadd -f histograms_tau.root histograms.root ../gamma3/histograms.root
cd ..
python3 plotter_tau.py ${masses[@]}


echo All done ":)"