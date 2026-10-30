# Rzut ukośny bez oporu i z oporem powietrza

*Raport obliczeniowy. Porównanie rozwiązania analitycznego z numerycznym modelem oporu aerodynamicznego.*

**Parametry początkowe:** $v_0=37\,\mathrm{m/s}$, $\alpha=52^\circ$.  
**Metoda:** wyprowadzenie równań, obliczenia w SI, całkowanie numeryczne metodą Rungego-Kutty 4. rzędu.  
**Charakter wyników:** symulacja modelowa, nie pomiar laboratoryjny.

Rzut ukośny jest ruchem w dwóch kierunkach jednocześnie. W modelu idealnym pozioma składowa prędkości pozostaje stała, a pionowa zmienia się pod wpływem grawitacji. Opór powietrza sprzęga obie składowe i zmienia tor z paraboli w krzywą niesymetryczną.

## Spis treści

1. [Cel, założenia i dane](#1-cel-założenia-i-dane)
2. [Rzut bez oporu powietrza](#2-rzut-bez-oporu-powietrza)
3. [Rzut z oporem powietrza](#3-rzut-z-oporem-powietrza)
4. [Wyniki i wykresy](#4-wyniki-i-wykresy)
5. [Kod i odtwarzanie obliczeń](#5-kod-i-odtwarzanie-obliczeń)
6. [Wnioski i ograniczenia](#6-wnioski-i-ograniczenia)

---

## 1. Cel, założenia i dane

Celem jest wyznaczenie toru, czasu lotu, zasięgu, wysokości maksymalnej oraz zmian prędkości i energii. W pierwszej części pomijamy opór, co pozwala otrzymać wzory analityczne. W drugiej uwzględniamy opór kwadratowy i rozwiązujemy równania numerycznie.

Początek układu współrzędnych znajduje się w miejscu wyrzutu. Oś $x$ jest pozioma, oś $y$ skierowana pionowo w górę, a poziom gruntu ma równanie $y=0$. Kąt $\alpha$ mierzymy od dodatniego zwrotu osi $x$.

### Założenia modelu

- Przyspieszenie grawitacyjne jest stałe: $g=9{,}81\,\mathrm{m/s^2}$.
- Wyrzut i lądowanie odbywają się na tej samej wysokości.
- Nie ma wiatru, obrotu piłki ani siły Magnusa.
- Model z oporem opisuje piłkę jako kulę o stałym współczynniku oporu; gęstość powietrza nie zmienia się z wysokością.
- Opór jest proporcjonalny do kwadratu szybkości. Nie jest to model oporu liniowego dla bardzo małych prędkości.

| Symbol | Znaczenie | Wartość | Jednostka SI |
| :--- | :--- | ---: | :---: |
| $v_0$ | szybkość początkowa | $37$ | $\mathrm{m/s}$ |
| $\alpha$ | kąt wyrzutu | $52^\circ$ ($0{,}9076\,\mathrm{rad}$) | radian (kąt) |
| $g$ | przyspieszenie grawitacyjne | $9{,}81$ | $\mathrm{m/s^2}$ |
| $m$ | masa modelowanej piłki | $0{,}145$ | $\mathrm{kg}$ |
| $\rho$ | gęstość powietrza | $1{,}225$ | $\mathrm{kg/m^3}$ |
| $C_D$ | bezwymiarowy współczynnik oporu | $0{,}47$ | $1$ |
| $d$ | średnica kuli | $0{,}073$ | $\mathrm{m}$ |
| $A$ | pole przekroju kuli | $0{,}004185$ | $\mathrm{m^2}$ |
| $\Delta t$ | krok całkowania numerycznego | $0{,}002$ | $\mathrm{s}$ |

Wartości $m$, $d$ i $C_D$ są przykładowymi parametrami piłki zbliżonej rozmiarem do baseballowej. Zmiana tych danych zmieni wyniki z oporem, ale nie zmieni rozwiązania idealnego.

## 2. Rzut bez oporu powietrza

### 2.1. Rozkład prędkości początkowej

Prędkość początkową rozkładamy na składowe poziomą i pionową:

$$
v_{x0}=v_0\cos\alpha,\qquad v_{y0}=v_0\sin\alpha.
$$

Dla podanych danych:

$$
v_{x0}=37\cos 52^\circ=22{,}779\,\mathrm{m/s},\qquad
v_{y0}=37\sin 52^\circ=29{,}156\,\mathrm{m/s}.
$$

### 2.2. Równania ruchu i ich rozwiązanie

Bez oporu jedyną siłą jest ciężar $\vec F_g=(0,-mg)$. Z drugiej zasady Newtona, $m\vec a=\vec F$, dostajemy:

$$
\ddot{x}=0,\qquad \ddot{y}=-g.
$$

Całkujemy równania, wykorzystując warunki początkowe $x(0)=y(0)=0$, $\dot{x}(0)=v_{x0}$ i $\dot{y}(0)=v_{y0}$:

$$
x(t)=v_{x0}t=v_0\cos\alpha\,t,
\qquad
y(t)=v_{y0}t-\frac{gt^2}{2}=v_0\sin\alpha\,t-\frac{gt^2}{2}.
$$

Różniczkując po czasie, otrzymujemy składowe prędkości:

$$
v_x(t)=v_{x0},\qquad v_y(t)=v_{y0}-gt.
$$

Tor otrzymujemy przez wyeliminowanie czasu: z pierwszego równania $t=x/v_{x0}$. Po podstawieniu do $y(t)$:

$$
y(x)=x\tan\alpha-\frac{g x^2}{2v_0^2\cos^2\alpha}.
$$

Jest to parabola. Jej ramiona są skierowane w dół, ponieważ współczynnik przy $x^2$ jest ujemny.

### 2.3. Czas lotu, zasięg i wysokość

W najwyższym punkcie $v_y=0$, zatem:

$$
t_{\max}=\frac{v_{y0}}{g}=2{,}972\,\mathrm{s},\qquad
H= y(t_{\max})=\frac{v_{y0}^2}{2g}=43{,}328\,\mathrm{m}.
$$

Przy lądowaniu $y(t_f)=0$. Rozwiązania to $t=0$ (chwila wyrzutu) oraz $t_f=2v_{y0}/g$:

$$
t_f=\frac{2v_0\sin\alpha}{g}=5{,}944\,\mathrm{s}.
$$

Zasięg poziomy wynosi $R=x(t_f)$:

$$
R=v_{x0}t_f=\frac{v_0^2\sin(2\alpha)}{g}=135{,}406\,\mathrm{m}.
$$

Wzór na zasięg wykorzystuje tożsamość $2\sin\alpha\cos\alpha=\sin(2\alpha)$. Dla stałej prędkości początkowej i jednakowej wysokości startu oraz lądowania maksimum zasięgu bez oporu wypada przy $\alpha=45^\circ$. Nie należy przenosić tego wniosku bez zmian na model z oporem.

## 3. Rzut z oporem powietrza

### 3.1. Siła oporu i równania

Dla oporu aerodynamicznego proporcjonalnego do kwadratu szybkości przyjmujemy:

$$
\vec F_D=-\frac{1}{2}\rho C_D A\lvert\vec v\rvert\vec v,
\qquad
A=\pi\left(\frac{d}{2}\right)^2,
\qquad
\lvert\vec v\rvert=\sqrt{v_x^2+v_y^2}.
$$

Znak minus oznacza, że siła jest skierowana przeciwnie do prędkości. Składowe drugiej zasady Newtona mają postać:

$$
\ddot{x}=-k\sqrt{\dot{x}^2+\dot{y}^2}\,\dot{x},\qquad
\ddot{y}=-g-k\sqrt{\dot{x}^2+\dot{y}^2}\,\dot{y},
\qquad
k=\frac{\rho C_D A}{2m}.
$$

Dla przyjętych parametrów $k\approx0{,}00831\,\mathrm{m^{-1}}$. Równania są sprzężone: pozioma składowa oporu zależy także od $v_y$, a pionowa także od $v_x$. Nie ma tu prostego rozwiązania parabolicznego.

W zapisie pierwszego rzędu wprowadzamy stan $\mathbf{s}=(x,y,v_x,v_y)$:

$$
\dot{\mathbf{s}}=\left(v_x,\ v_y,
-k\lvert\vec v\rvert v_x,
-g-k\lvert\vec v\rvert v_y\right).
$$

### 3.2. Całkowanie numeryczne

Równanie stanu rozwiązujemy metodą Rungego-Kutty 4. rzędu (RK4). Dla $\dot{\mathbf{s}}=\mathbf{f}(t,\mathbf{s})$ i kroku $h=\Delta t$ pojedynczy krok ma postać:

$$
\begin{aligned}
\mathbf{k}_1&=\mathbf{f}(t_n,\mathbf{s}_n),\\
\mathbf{k}_2&=\mathbf{f}(t_n+h/2,\mathbf{s}_n+h\mathbf{k}_1/2),\\
\mathbf{k}_3&=\mathbf{f}(t_n+h/2,\mathbf{s}_n+h\mathbf{k}_2/2),\\
\mathbf{k}_4&=\mathbf{f}(t_n+h,\mathbf{s}_n+h\mathbf{k}_3),\\
\mathbf{s}_{n+1}&=\mathbf{s}_n+\frac{h}{6}(\mathbf{k}_1+2\mathbf{k}_2+2\mathbf{k}_3+\mathbf{k}_4).
\end{aligned}
$$

Obliczenia kończymy po przekroczeniu poziomu $y=0$ w dół. Położenie chwili lądowania interpolujemy liniowo między ostatnim punktem nad ziemią i pierwszym pod ziemią. Krok $0{,}002\,\mathrm{s}$ jest na tyle mały, by dać stabilne wartości do przedstawionego porównania.

### 3.3. Bilans energii

Energia mechaniczna to suma energii kinetycznej i potencjalnej grawitacji:

$$
E_m=E_k+E_p=\frac{m}{2}(v_x^2+v_y^2)+mgy.
$$

W modelu bez oporu jej pochodna po czasie jest równa zeru. Z oporem:

$$
\frac{dE_m}{dt}=\vec F_D\cdot\vec v
=-\frac{1}{2}\rho C_D A\lvert\vec v\rvert^3\leq 0.
$$

Ujemna pochodna oznacza, że część energii mechanicznej jest rozpraszana, między innymi w postaci ciepła. Dlatego krzywa energii mechanicznej w modelu z oporem opada.

## 4. Wyniki i wykresy

| Wielkość | Bez oporu | Z oporem kwadratowym |
| :--- | ---: | ---: |
| Czas lotu $t_f$ | $5{,}944\,\mathrm{s}$ | $4{,}917\,\mathrm{s}$ |
| Zasięg $R$ | $135{,}406\,\mathrm{m}$ | $73{,}582\,\mathrm{m}$ |
| Wysokość maksymalna $H$ | $43{,}328\,\mathrm{m}$ | $29{,}853\,\mathrm{m}$ |
| Składowa pozioma prędkości na starcie | $22{,}779\,\mathrm{m/s}$ | $22{,}779\,\mathrm{m/s}$ |
| Składowa pionowa prędkości na starcie | $29{,}156\,\mathrm{m/s}$ | $29{,}156\,\mathrm{m/s}$ |

W tym przykładzie opór skraca zasięg o około $45{,}7\%$ i wysokość maksymalną o około $31{,}1\%$ względem modelu idealnego. Wyniki z oporem zależą od masy, rozmiaru, współczynnika oporu i gęstości ośrodka.

![Wykres porównujący paraboliczny tor bez oporu oraz krótszy i niesymetryczny tor z oporem aerodynamicznym](files/rzut_ukosny/tor_lotu.png)

*Rysunek 1. Tor lotu. Wysokości i zasięgi odpowiadają wartościom w tabeli; wykresy są poglądowymi wizualizacjami obliczeń modelowych.*

![Wykres położenia poziomego i wysokości w czasie, porównujący oba modele](files/rzut_ukosny/polozenie_czas.png)

*Rysunek 2. W modelu idealnym $x(t)$ rośnie liniowo. Z oporem pozioma prędkość maleje; wysokość w obu modelach najpierw rośnie, a następnie spada.*

![Wykres składowych prędkości i szybkości w czasie dla rzutu bez oporu oraz z oporem](files/rzut_ukosny/predkosc_czas.png)

*Rysunek 3. Bez oporu $v_x$ jest stałe, a $v_y$ zmienia się liniowo. Opór zmniejsza obie składowe, a szybkość przestaje być symetryczna względem wierzchołka toru.*

![Wykres energii kinetycznej, potencjalnej i całkowitej z zachowaniem energii w modelu idealnym i spadkiem energii przy oporze](files/rzut_ukosny/energia.png)

*Rysunek 4. Bez oporu energia mechaniczna jest stała. Opór wykonuje ujemną pracę i zmniejsza energię mechaniczną układu.*

## 5. Kod i odtwarzanie obliczeń

Pełny generator wykresów znajduje się w [generuj_wykresy.py](generuj_wykresy.py). Korzysta z NumPy i Matplotlib, wykonuje całkowanie RK4, zapisuje wykresy w formacie PNG oraz sprawdza podstawowe warunki modelu.

Najważniejszy fragment modelu oporu i całkowania:

```python
def acceleration(state):
    vx, vy = state[2], state[3]
    speed = np.hypot(vx, vy)
    drag_per_mass = 0.5 * RHO * CD * AREA * speed / MASS
    return np.array([
        vx,
        vy,
        -drag_per_mass * vx,
        -G - drag_per_mass * vy,
    ])

k1 = acceleration(state)
k2 = acceleration(state + 0.5 * DT * k1)
k3 = acceleration(state + 0.5 * DT * k2)
k4 = acceleration(state + DT * k3)
next_state = state + DT * (k1 + 2*k2 + 2*k3 + k4) / 6
```

Generator wymaga Pythona oraz pakietów `numpy` i `matplotlib`. Po ich zainstalowaniu uruchom w katalogu raportu:

```powershell
python -m pip install numpy matplotlib
python generuj_wykresy.py
```

Wykresy zostaną zapisane w `files/rzut_ukosny/`. Zmieniając stałe na początku skryptu, można sprawdzić wpływ masy i parametrów aerodynamicznych; aby zachować zgodność raportu, $v_0$ i $\alpha$ należy pozostawić równe wartościom wynikającym z indeksu.

## 6. Wnioski i ograniczenia

1. Bez oporu ruch poziomy jest jednostajny, a pionowy jednostajnie przyspieszony. Po eliminacji czasu tor ma kształt paraboli.
2. W modelu oporu kwadratowego siła zależy od wektora prędkości, więc równania poziome i pionowe są sprzężone i wymagają całkowania numerycznego.
3. Dla przyjętego modelu piłki opór wyraźnie ogranicza zasięg i wysokość oraz powoduje spadek energii mechanicznej.
4. Wyniki nie opisują dokładnie rzeczywistej piłki: pomijają wiatr, rotację, zmianę gęstości powietrza oraz zależność $C_D$ od prędkości i orientacji. Są porównaniem dwóch jasno określonych modeli.

## Dalsza lektura

- [OpenStax: Projectile Motion](https://openstax.org/books/college-physics-2e/pages/3-4-projectile-motion)
- [NASA Glenn: The Drag Equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/)