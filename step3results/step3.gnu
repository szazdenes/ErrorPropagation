reset

set terminal postscript enhanced portrait font 'Times-Roman,12'
set output 'step3_north_new.ps'

set style line 1 lc rgb 'grey' lt 1 pt 7 lw 1
set style line 2 lc rgb 'red' lt 1 pt 7 lw 2
set style line 3 lc rgb 'black' lt 1 pt 7 lw 1

filename1 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0621_am_best.csv"
filename2 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0621_am_worst.csv"
filename3 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0621_pm_best.csv"
filename4 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0621_pm_worst.csv"
filename5 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0922_am_best.csv"
filename6 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0922_am_worst.csv"
filename7 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0922_pm_best.csv"
filename8 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0922_pm_worst.csv"
filename9 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0621_am_secondworst.csv"
filename10 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0621_pm_secondworst.csv"
filename11 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0922_am_secondworst.csv"
filename12 = "/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error_0922_pm_secondworst.csv"

set multiplot layout 4,3

#################################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename1 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename9 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename2 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename3 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename10 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename4 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename5 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename11 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename6 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename7 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename12 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

set xrange [0:60]
set yrange [-65:65]
set xlabel 'Sun elevation ({^o})'
set ylabel 'North error ({/Symbol w}_{North}) ({^o})
set zeroaxis
set tics front
plot filename8 u 1:2:3 ls 1 w filledcurve notitle, '' u 1:2 ls 3 w l notitle, '' u 1:3 ls 3 w l notitle

###########################################

unset multiplot

