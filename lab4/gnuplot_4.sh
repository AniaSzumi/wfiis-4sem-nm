set terminal post enhanced colour solid font 25  # wybor formatu, w jakim chcemy utworzyc wynikowy rysunek
set encoding iso_8859_2

set style line 1 pt 7 lc rgb "red"           
set style line 2 pt 7 lc rgb "black" 
set style line 3 pt 7 lc rgb "dark-spring-green"               
set style line 4 pt 7 lc rgb "goldenrod" 
set style line 5 pt 7 lc rgb "turquoise"               
set style line 6 pt 7 lc rgb "royalblue" 

# eval.dat -- plik z pierwiastkami wartosci wlasnych. Zawiera 7 kolumn: parametr alfa oraz 6 odpowiadajacych mu pierwiastkow wartosci wlasnych.
set output "eval.eps" # Pierwiastki wartosci wlasnych w funkcji parametru alpha
set xlabel "{/Symbol a}"
set ylabel "{/Symbol w}({/Symbol a})"
set key outside horizontal top center
set multiplot
plot for [k=0:5] "eval.txt" u 1:(column(k+2)) w lp ls k+1 t "{/Symbol w}_".k.""
# Wstawiony rysunek wewnetrzny:
set origin 0.25,0.2
set size 0.7,0.7
unset key
set ylabel offset 1.5,0
set logscale x
set xtics ("2" 2, "10" 10, "100" 100)
replot
unset multiplot

reset
set style line 1 pt 7 lw 4 lc rgb "red"           
set style line 2 pt 7 lw 4 lc rgb "black" 
set style line 3 pt 7 lw 4 lc rgb "dark-spring-green"               
set style line 4 pt 7 lw 4 lc rgb "goldenrod" 
set style line 5 pt 7 lw 4 lc rgb "turquoise"               
set style line 6 pt 7 lw 4 lc rgb "royalblue" 

# evec.dat -- plik z wektorami wlasnymi. Zawiera 7 kolumn: wspolrzedna x_i oraz 6 odpowiadajacych mu elementow kolejnych wektorow. 
# Pierwsza seria danych: dla alpha = 0, druga: dla alpha = 100.
# Serie danych musza byc oddzielone dwiema pustymi liniami.
set output "evec0.eps" # Wektory wlasne w funkcji x dla alpha = 0
set xlabel "x"
set ylabel "u(x)"
set key outside horizontal top center
nrSerii=0
set label 1 "{/Symbol a} = 0" at screen 0.1, 0.93 center tc rgb "brown"
plot [-5:5][-.1157:.1157] for [k=0:5] "evec.txt" i nrSerii u 1:(column(k+2)) w l ls k+1 t "u_".k."(x)"

set output "evec100.eps" # Wektory wlasne w funkcji x dla alpha = 100
nrSerii=1
set label 1 "{/Symbol a} = 100" at screen 0.1, 0.93 center tc rgb "brown"
plot [-5:5][-.1157:.1157] for [k=0:5] "evec.txt" i nrSerii u 1:(column(k+2)) w l ls k+1 t "u_".k."(x)"
