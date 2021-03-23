set terminal post enhanced colour solid font 15 # wybor formatu, w jakim chcemy utworzyc wynikowy rysunek
set encoding iso_8859_2  # kodowanie

# Ustawienie styli linii (style 1 i 2):
set style line 1 pt 7 ps 0.6 lc rgb "red"               
set style line 2 pt 7 ps 0.6 lc rgb "blue" 


set output "rys.eps"      # Nazwa obrazka
set size 1, 0.65         # Proporcje wymiarow obrazka
set key right center      # Polozenie legendy
set xlabel "k (numer iteracji)"                      # Podpis osi OX
set ylabel "||r_k||_2" tc rgb "red"       # Podpis LEWEJ osi OY
set y2label "||x_k||_2" tc rgb "blue" rotate by -90  # Podpis PRAWEJ osi OY
set ytics nomirror        # Uniezaleznienie znacznikow skali osi OY od drugiej osi OY
set y2tics nomirror       # Analogicznie.
set logscale y            # Skala logarytmiczna na LEWEJ osi OY
set ytics tc rgb "red"    # Kolor podpisow znacznikow LEWEJ osi OY
set y2tics tc rgb "blue"  # Kolor podpisow znacznikow PRAWEJ osi OY
set grid                  # Wlaczenie widocznosci siatki
set y2range [4695:]       # Ustawienie dolnej granicy przedzialu PRAWEJ osi OY

# Komenda rysowania dwoch serii danych z pliku tekstowego:
plot "out.txt" u 1:2 w lp ls 1 t "||r_k||_2,", \
"" u 1:5 w lp ls 2 t "||x_k||_2" axes x1y2

# INFORMACJE NA TEMAT PLIKU Z DANYMI:
#
# Plik: out.dat zawiera wyniki
#          - W pierwszej kolumnie: k        -- numer iteracji
#          - W drugiej kolumnie:   ||r_k||  -- norma z r_k
#          - W trzeciej kolumnie:  alfa_k
#          - W czwartej kolumnie:  beta_k
#          - W piatej kolumnie:  ||x_k||  -- norma z x_k

# WYJASNIENIE OPCJI KOMENDY plot:
#          u 1:2 == using 1:2 -- 1. kolumna danych zostanie uznana za wartosci osi OX, 2. kolumna: OY (u 1:5 -- analogicznie przy uzyciu 1. i 5. kolumny danych)
#          ls 1  == linestyle 1 -- wybor stylu linii
#          t "tekst" == title "tekst" -- "tekst" bedzie podpisem na legendzie
#          axes x1y2 -- uzycie PRAWEJ osi OY (bez tej opcji Gnuplot domyslnie wybierze LEWA os OY).
