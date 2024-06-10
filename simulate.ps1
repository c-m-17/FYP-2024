# simulate.ps1

$pythonScriptPath = "C:\Users\cerys\OneDrive\Documents\CivEng\!FYP\FYP-2024\main.py"

for($d=2; $d -le 12; $d++){

    $diam = $d * 0.5

    for($g=1; $g -le 4; $g++){

        $geoms = "inputs\geometries-0$g.csv"

        for($l=1; $l -le 4; $l++){

            $loads = "inputs\loading-0$l.csv"

            Start-Process -FilePath $pythonScriptPath -ArgumentList "-g $geoms -l $loads -d $diam" -Wait
        }
    }
}