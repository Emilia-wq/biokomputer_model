# Biocomputer DNA/RNA - Molekularna Bramka Logiczna

## Przegląd projektu

Projekt implementuje molekularny biokomputer wykorzystujący bramki logiczne RNA (AND, NOT) do detekcji komórek nowotworowych poprzez specyficzne markery mikroRNA. System generuje sygnał wyjściowy wyłącznie przy precyzyjnej kombinacji biomarkerów miRNA, umożliwiając rozróżnienie komórek nowotworowych od zdrowych.

**Funkcja logiczna:** `Signal = (miR-21 AND miR-155) AND NOT miR-34`

## Struktura projektu
```
├── NUPACK_Generator.ipynb      # Python: Analiza NUPACK i generacja SBML
├── BioTest.m                   # MATLAB: Symulacje SimBiology i analiza
├── biocomputer_tube_design.xml # Wygenerowany model SBML
├── Raport.pdf                  # Szczegółowy raport techniczny
└── README.md
```

## Metodologia

### Faza 1: Projektowanie sekwencji (Python + NUPACK)

- **Narzędzie:** NUPACK 4.0 - analiza termodynamiczna
- **Cel:** Zaprojektowanie optymalnych sekwencji bramek RNA z minimalną strukturą drugorzędową (folding) i crosstalkiem
- **Metoda:** Algorytm multi-tube design z ograniczeniami twardymi i miękkimi
- **Wynik:** Model SBML z termodynamicznie zwalidowanymi parametrami kinetycznymi

**Kluczowe funkcjonalności:**
- Projektowanie przełączników toehold dla mechanizmu inhibicji miR-34
- Korekta folding penalty (redukcja k_on na podstawie Delta G)
- Kompleksowa analiza crosstalk (6 par off-target)
- Obliczanie współczynnika ortogonalności

### Faza 2: Symulacje dynamiczne (MATLAB + SimBiology)

- **Narzędzie:** MATLAB R2023a z SimBiology Toolbox
- **Solver:** Deterministyczny (ode15s) dla stabilności numerycznej
- **Zakres analiz:**
  - Weryfikacja tabeli prawdy (8 kombinacji wejść)
  - Kwantyfikacja stosunku ON/OFF
  - EC50 (wrażliwość aktywacji)
  - IC50 (efektywność inhibicji)
  - Testy odporności (przeciek, crosstalk)

## Tło biologiczne

### Biomarkery mikroRNA

| miRNA | Funkcja | Ekspresja w raku |
|-------|---------|------------------|
| **miR-21** | Onkogen | Nadekspresja |
| **miR-155** | Onkogen | Nadekspresja |
| **miR-34a** | Supresor | Obniżona/wyciszona |

**Logika diagnostyczna:**
- `(miR-21 HIGH AND miR-155 HIGH AND miR-34 LOW)` oznacza **Sygnatura nowotworowa wykryta**
- Pozostałe kombinacje oznaczają **Profil prawidłowy/niejednoznaczny**

## Instalacja i użytkowanie

### Wymagania

**Środowisko Python:**
```bash
pip install nupack python-libsbml numpy
```

**MATLAB:**
- MATLAB R2023a lub nowszy
- SimBiology Toolbox

## Opis plików

### BioTest.ipynb

Notebook Python implementujący:
- Klasę `AdvancedNUPACKAnalyzer`
  - `design_gates_with_tube_design()`: Optymalizacja sekwencji multi-tube
  - `analyze_secondary_structure()`: Predykcja struktur MFE + folding penalties
  - `analyze_binding()`: Obliczanie K_d termodynamicznego przez tube_analysis
  - `analyze_crosstalk()`: Kwantyfikacja wiązań off-target
- Funkcje eksportu SBML z termodynamicznie skorygowanymi parametrami

### BioTest.m

Skrypt MATLAB realizujący:
1. Import modelu SBML
2. Weryfikacja tabeli prawdy (8 scenariuszy)
3. Analiza odpowiedzi dynamicznej (stan ON)
4. Krzywe dawka-odpowiedź (EC50, IC50)
5. Testy stress odporności (przeciek, crosstalk)
6. Wizualizacja i eksport danych

### biocomputer_tube_design.xml

Model SBML Level 3 zawierający:
- 13 gatunków chemicznych (miRNA, bramki, intermediaty, sygnał)
- 23 reakcje (wiązanie, bramki AND/NOT, degradacja)
- Parametry kinetyczne z NUPACK (k_on, k_off)

## Ograniczenia i przyszłe kierunki

### Obecne ograniczenia

1. **Tylko symulacje deterministyczne** - Efekty stochastyczne (fluktuacje małej liczby kopii) nieuwzględnione
2. **Uproszczona degradacja** - Jednolite k_deg = 10⁻³ s⁻¹ dla wszystkich gatunków
3. **Pojedynczy kompartment** - Brak niejednorodności przestrzennej i ograniczeń dyfuzyjnych
4. **Umiarkowany crosstalk** - miR-21/Gate155 (K_d = 108 nM) może obniżać specyficzność przy wysokich stężeniach

### Proponowane usprawnienia

- **Symulacje stochastyczne** (algorytm SSA Gillespiego) dla gatunków o niskiej abundancji
- **Walidacja eksperymentalna** z użyciem in vitro transkrybowanego RNA i odczytu fluorescencyjnego
- **Kaskady wielobramkowe** dla logiki wyższego rzędu (3+ wejścia)
- **Przeprojektowanie sekwencji** w celu eliminacji umiarkowanego crosstalk (cel: wszystkie K_d > 1 µM)
