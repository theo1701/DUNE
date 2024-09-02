set table "conts1"
set view map
set samples 25,25
set isosamples 26,26

unset surface
set contour base

set cntrparam bspline
set cntrparam levels discrete 2.30

pi = acos(-1.0)
plotdcp(x) = x*180/pi 
plott13(x) = sin(2.0*x)*sin(2.0*x)

splot "dchsq_vs_dcptest+th23test_NHtrue_NOvA+T2Kbestfit_NHtest_DUNE5+5.dat" u 1:2:3 title 'NH test/1 sigma' #ls 1

unset table
###########################################

set table "conts2"
set view map
set samples 25,25
set isosamples 26,26

unset surface
set contour base

set cntrparam bspline
set cntrparam levels discrete 2.30

pi = acos(-1.0)
plotdcp(x) = x*180/pi 
plott13(x) = sin(2.0*x)*sin(2.0*x)

splot "dchsq_vs_dcptest+th23test_NHtrue_NOvA+T2Kbestfit_IHtest_DUNE5+5.dat" u 1:2:3 title 'IH test/1 sigma' #ls 1

unset table
###########################################

set table "conts3"
set view map
set samples 25,25
set isosamples 26,26

unset surface
set contour base

set cntrparam bspline
set cntrparam levels discrete  2.30

pi = acos(-1.0)
plotdcp(x) = x*180/pi 
plott13(x) = sin(2.0*x)*sin(2.0*x)

splot "dchsq_vs_dcptest+th23test_IHtrue_NOvA+T2Kbestfit_NHtest_DUNE5+5.dat" u 1:2:3 title 'NH test/1 sigma' #ls 1

unset table
###########################################


set table "conts4"
set view map
set samples 25,25
set isosamples 26,26

unset surface
set contour base

set cntrparam bspline
set cntrparam levels discrete  2.30

pi = acos(-1.0)
plotdcp(x) = x*180/pi 
plott13(x) = sin(2.0*x)*sin(2.0*x)

splot "dchsq_vs_dcptest+th23test_IHtrue_NOvA+T2Kbestfit_IHtest_DUNE5+5.dat" u 1:2:3 title 'IH test/1 sigma' #ls 1

unset table
###########################################
set table "conts5"
set view map
set samples 25,25
set isosamples 26,26

unset surface
set contour base

set cntrparam bspline
set cntrparam levels discrete  11.83

pi = acos(-1.0)
plotdcp(x) = x*180/pi 
plott13(x) = sin(2.0*x)*sin(2.0*x)

splot "dchsq_vs_dcptest+th23test_NHtrue_NOvA+T2Kbestfit_NHtest_DUNE5+5.dat" u 1:2:3 title 'NH test/3 sigma' #ls 1

unset table
###########################################


set table "conts6"
set view map
set samples 25,25
set isosamples 26,26

unset surface
set contour base

set cntrparam bspline
set cntrparam levels discrete  11.83

pi = acos(-1.0)
plotdcp(x) = x*180/pi 
plott13(x) = sin(2.0*x)*sin(2.0*x)

splot "dchsq_vs_dcptest+th23test_NHtrue_NOvA+T2Kbestfit_IHtest_DUNE5+5.dat" u 1:2:3 title 'IH test/3 sigma' #ls 1

unset table
###########################################
set table "conts7"
set view map
set samples 25,25
set isosamples 26,26

unset surface
set contour base

set cntrparam bspline
set cntrparam levels discrete  11.83

pi = acos(-1.0)
plotdcp(x) = x*180/pi 
plott13(x) = sin(2.0*x)*sin(2.0*x)

splot "dchsq_vs_dcptest+th23test_IHtrue_NOvA+T2Kbestfit_NHtest_DUNE5+5.dat" u 1:2:3 title 'NH test/3 sigma' #ls 1

unset table
###########################################


set table "conts8"
set view map
set samples 25,25
set isosamples 26,26

unset surface
set contour base

set cntrparam bspline
set cntrparam levels discrete  11.83

pi = acos(-1.0)
plotdcp(x) = x*180/pi 
plott13(x) = sin(2.0*x)*sin(2.0*x)

splot "dchsq_vs_dcptest+th23test_IHtrue_NOvA+T2Kbestfit_IHtest_DUNE5+5.dat" u 1:2:3 title 'IH test/3 sigma' #ls 1

unset table
####################################################################################################################
#set xlabel "{/Symbol d}_{CP} (true)" 
#set xrange [ 0:1.2 ] noreverse nowriteback
#set logscale y
#set ylabel "Probability"
#set yrange [0 : 25] noreverse nowriteback
#set zrange [ -1 : 50 ] noreverse nowriteback
#set xtics 0,0.2,1.2
#set ytics 0,5,25
set grid xtics 
set grid ytics
#set key at 170, 0.9*10**7

AB=0.40
CD=0.2
EF=0.1

P1=0.12 ; Q1= EF
P2=0.12 ; Q2= EF+AB
P3=0.12+AB ; Q3= EF+AB
P4=0.12+AB; Q4= EF


set bmargin CD
set tmargin CD
set lmargin CD
set rmargin CD


set style line 1 dashtype 1 lc rgb "black" lw 1
set style line 2 dashtype 2 lc rgb "black" lw 1
set style line 3 dashtype 4 lc rgb "black" lw 1
set style line 4 dashtype 8 lc rgb "black" lw 1
set style line 5 dashtype 1 lc rgb "#228B22" lw 1
set style line 6 dashtype 2 lc rgb "#228B22" lw 1
set style line 7 dashtype 4 lc rgb "#228B22" lw 1
set style line 8 dashtype 8 lc rgb "#228B22" lw 1
set style line 9 dashtype 1 lc rgb "red" lw 1
set style line 10 dashtype 2 lc rgb "red" lw 1
#set style line 11 dashtype  4 lc rgb "red" lw 1
set style line 12 dashtype 8 lc rgb "red" lw 1
set style line 13 dashtype 1 lc rgb "blue" lw 1
set style line 14 dashtype 2 lc rgb "blue" lw 1
#set style line 15 dashtype  4 lc rgb "blue" lw 1
set style line 16 dashtype 8 lc rgb "blue" lw 1

set style line 1 lc rgb 'black' pt 5   # square
set style line 2 lc rgb 'black' pt 7   # circle
set style line 3 lc rgb 'black' pt 9   # triangle
set style increment user 
###############################################################

set terminal postscript color enhanced 
set term postscript eps font "Times-Roman,18"

set output "dchsq_vs_dcptest+th23test_NHtrue_IHtrue_NOvA+T2Kbestfit_NHtest_IHtest_DUNE5+5.eps"


#######################################################################################################
#set size square 1
#set origin 0,0
set multiplot layout 1,2 columnsfirst scale 1,2



#######################################################################################################
set size AB,AB
set origin P2,Q2

set label "NH true, NH test" at -80, 0.65
set yrange [ 0.3: 0.7] noreverse nowriteback

set xrange [-180 : 180] noreverse nowriteback
#set zlabel "Z axis" 
#set zlabel  offset character 1, 0, 0 font "" textcolor lt -1 norotate
#set zrange [ -1 : 50 ] noreverse nowriteback

set ytics (""0.3, "0.4"0.4, "0.5"0.5, "0.6"0.6, "0.7"0.7)
set xtics (""-180, ""-90, ""0, ""90, ""180)
#set grid xtics 
#set grid ytics
set mxtics 5
set mytics 5
#set label "{/Symbol d}_{cp}(test)" at -130, 0.25
set key at 175,0.48



plot "conts1" w l ls 9 t "1 {/Symbol s}",  "conts5" w l ls 10 t "3 {/Symbol s}", "bestfit_NH.dat" w p ls 2 t "LBL NH best fit (-170^{o}, 0.58)"
unset label
unset xlabel
unset ylabel
unset xrange
unset yrange
unset key
unset title

#######################################################################################################
set size AB,AB
set origin P3,Q3

set label "NH true, IH test" at -80, 0.65
set yrange [ 0.3: 0.7] noreverse nowriteback

set xrange [-180 : 180] noreverse nowriteback
#set zlabel "Z axis" 
#set zlabel  offset character 1, 0, 0 font "" textcolor lt -1 norotate
#set zrange [ -1 : 50 ] noreverse nowriteback

set ytics (""0.3, ""0.4, ""0.5, ""0.6, ""0.7)
set xtics (""-180, ""-90, ""0, ""90, ""180)
#set grid xtics 
#set grid ytics
set mxtics 5
set mytics 5
#set label "{/Symbol d}_{cp}(test)" at -130, 0.25
set key at 210,0.48



plot "conts2" w l ls 13 t "1 {/Symbol s}",  "conts6" w l ls 14 t "3 {/Symbol s}",  "bestfit_NOvA+T2K_nu+anu_app+disapp.dat" w p ls 3 title "LBL best fit (-90^{o}, 0.58)"

unset label
unset xlabel
unset ylabel
unset xrange
unset yrange
unset key
unset title
#######################################################################################################
#######################################################################################################
set size AB,AB
set origin P1,Q1

set label "IH true, NH test" at -80, 0.65
set yrange [ 0.3: 0.7] noreverse nowriteback

set xrange [-180 : 180] noreverse nowriteback
#set zlabel "Z axis" 
#set zlabel  offset character 1, 0, 0 font "" textcolor lt -1 norotate
#set zrange [ -1 : 50 ] noreverse nowriteback

set ytics ("0.3"0.3, "0.4"0.4, "0.5"0.5, "0.6"0.6, "0.7/0.3"0.7)
set xtics ("-180"-180, "-90"-90, "0"0, "90"90, "180/-180"180)
#set grid xtics 
#set grid ytics
set mxtics 5
set mytics 5
#set label "{/Symbol d}_{cp}" at -130, 0.25
set key at 210,0.48

set label "{/Symbol d}_{CP}" at 175,0.22

set label "sin^{2}{/Symbol q}_{23}" at -290, 0.75
plot "conts3" w l ls 9 t "1 {/Symbol s}",  "conts7" w l ls 10 t "3 {/Symbol s}", "bestfit_NOvA+T2K_NHtest_nu+anu_app+disapp.dat" w p ls 1 title "LBL NH best fit (-170^{o}, 0.58)"
unset label
unset xlabel
unset ylabel
unset xrange
unset yrange
unset key
unset title

#######################################################################################################
set size AB,AB
set origin P4,Q4

set label "IH true, IH test" at -80, 0.65
set yrange [ 0.3: 0.7] noreverse nowriteback

set xrange [-180 : 180] noreverse nowriteback
#set zlabel "Z axis" 
#set zlabel  offset character 1, 0, 0 font "" textcolor lt -1 norotate
#set zrange [ -1 : 50 ] noreverse nowriteback

set ytics (""0.3, ""0.4, ""0.5, ""0.6, ""0.7)
set xtics (""-180, "-90"-90, "0"0, "90"90, "180"180)
#set grid xtics 
#set grid ytics
set mxtics 5
set mytics 5
#set label "{/Symbol d}_{cp}(test)" at -130, 0.25
set key at 179.99,0.48



plot "conts4" w l ls 13 t "1 {/Symbol s}",  "conts8" w l ls 14 t "3 {/Symbol s}", "bestfit_IH.dat" w p ls 2 t "LBL best fit (-90^{o}, 0.58)"
unset label
unset xlabel
unset ylabel
unset xrange
unset yrange
unset key
unset title
#######################################################################################################

unset multiplot

system "epstopdf dchsq_vs_dcptest+th23test_NHtrue_IHtrue_NOvA+T2Kbestfit_NHtest_IHtest_DUNE5+5.eps &"
system "gv dchsq_vs_dcptest+th23test_NHtrue_IHtrue_NOvA+T2Kbestfit_NHtest_IHtest_DUNE5+5.eps &"

unset table
###########################################


