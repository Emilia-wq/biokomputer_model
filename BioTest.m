%  Analiza molekularnego klasyfikatora DNA/RNA
%  Logika: Signal = (miR21 AND miR155) AND NOT miR34
%
%  Skrypt wczytuje model biokomputera (z pliku SBML) i przeprowadza
%  kompleksową analizę deterministyczną 
%  Autor: Emilia Romanowska - 3,4,6,8,9,10  Stanisław Wojtków - 1,2,5,7

clear; clc; close all;

%% 1. WCZYTANIE MODELU
sbml_filename = "biocomputer_tube_design.xml";

% Sprawdzenie, czy plik modelu istnieje
if ~isfile(sbml_filename)
    error("Nie znaleziono pliku SBML: %s", sbml_filename);
end

% Import modelu SBML do środowiska SimBiology
model = sbmlimport(sbml_filename);

%% 2. FUNKCJA POMOCNICZA DO SYMULACJI
function [sim_clean, max_signal] = simulate_scenario(model, miR21_conc, miR155_conc, miR34_conc, stop_time)

    %   Funkcja tworzy kopię modelu, ustawia podane stężenia 
    %   początkowe dla trzech wejść (miR21, miR155, miR34) i uruchamia 
    %   symulację za pomocą solwera 'ode15s'.
    %
    %   Wejścia:
    %     model (SimBiology.Model): Bazowy obiekt modelu.
    %     miR21_conc (double): Stężenie [M] wejścia miR21.
    %     miR155_conc (double): Stężenie [M] wejścia miR155.
    %     miR34_conc (double): Stężenie [M] wejścia miR34.
    %     stop_time (double): Czas końcowy symulacji [s].
    %
    %   Wyjścia:
    %     sim_clean (matrix): Oczyszczone dane symulacji (bez ujemnych wartości).
    %     max_signal (double): Maksymalna osiągnięta wartość 'Signal'.
    
    % Kopia obiektu modelu
    m = copyobj(model);
    
    % Ustawienie stężeń początkowych dla wejść
    s = m.Species;
    set(s(strcmp({s.Name}, 'miR21')), 'InitialAmount', miR21_conc);
    set(s(strcmp({s.Name}, 'miR155')), 'InitialAmount', miR155_conc);
    set(s(strcmp({s.Name}, 'miR34')), 'InitialAmount', miR34_conc);
    
    % Konfiguracja solwera deterministycznego
    cs = getconfigset(m, 'active');
    cs.SolverType = 'ode15s';
    cs.StopTime = stop_time;
    cs.SolverOptions.AbsoluteTolerance = 1e-20; % Tolerancja dla niskich stężeń
    cs.SolverOptions.RelativeTolerance = 1e-12;
    cs.SolverOptions.MaxStep = 0.5;
    
    % Uruchomienie symulacji z obsługą błędów
   try
        sim = sbiosimulate(m);
    catch ME 
        error('BŁĄD KRYTYCZNY w sbiosimulate: %s', ME.message);
    end
    
    % Oczyszczenie danych (stężenie nie może być ujemne)
    sim_clean = sim.Data;
    sim_clean(sim_clean < 0) = 0;
    
    % Walidacja wyjścia: Sprawdzenie, czy 'Signal' istnieje
    signal_idx = find(strcmp(sim.DataNames, 'Signal'), 1);
    
    if isempty(signal_idx)
        error('Nie znaleziono gatunku "Signal" w wynikach symulacji. Sprawdź, czy plik .xml jest poprawny.');
    end
    
    % Zwrócenie maksymalnej wartości sygnału
    max_signal = max(sim_clean(:, signal_idx));
end

%% 3. DEFINICJA TESTÓW I PARAMETRÓW GLOBALNYCH
% Definicja globalnych stałych
CONFIG.conc_high = 1e-6;  % 1 uM (Wysokie stężenie, stan '1')
CONFIG.conc_low = 0;      % 0 M  (Niskie stężenie, stan '0')
CONFIG.stop_time = 50;    % Domyślny czas symulacji [s]

% Definicja macierzy przypadków testowych
% Kolumna 4 definiuje oczekiwany stan (0=OFF, 1=ON, 2=INHIBITED)
% dla automatyzacji analizy w Sekcji 4.
test_cases = [
    % miR21  miR155  miR34  | Oczekiwany Stan (0=OFF, 1=ON, 2=INHIBITED)
      0      0       0        0
      0      0       1        0
      0      1       0        0
      0      1       1        0
      1      0       0        0
      1      0       1        0
      1      1       0        1  
      1      1       1        2  
];

% Weryfikacja tabeli prawdy
results = struct(); % Inicjalizacja struktury na wyniki
fprintf('Testowanie wszystkich %d kombinacji wejść\n', size(test_cases, 1));
fprintf('%-8s %-8s %-8s | %-15s | Oczekiwany Stan\n', 'miR21', 'miR155', 'miR34', 'Max Signal [M]');
fprintf('%s\n', repmat('-', 60, 1));

for i = 1:size(test_cases, 1)
    % Pobranie danych wejściowych z macierzy
    miR21_val = test_cases(i,1) * CONFIG.conc_high;
    miR155_val = test_cases(i,2) * CONFIG.conc_high;
    miR34_val = test_cases(i,3) * CONFIG.conc_high;
    
    % Pobranie oczekiwanego stanu 
    expected_state = test_cases(i, 4);
    
    % Uruchomienie symulacji
    [sim_data, max_sig] = simulate_scenario(model, miR21_val, miR155_val, miR34_val, CONFIG.stop_time);
    
    % Zapisanie wyników do struktury
    results(i).miR21 = test_cases(i,1);
    results(i).miR155 = test_cases(i,2);
    results(i).miR34 = test_cases(i,3);
    results(i).expected_state = expected_state;
    results(i).max_signal = max_sig;
    
    % Wyświetlanie statusu w konsoli
    status = 'OFF';
    if expected_state == 1
        status = 'ON';
    elseif expected_state == 2
        status = 'INHIBITED';
    end
    
    % Wyświetlenie wyników wiersza
    fprintf('%-8d %-8d %-8d | %-15.3e | %s\n', ...
        test_cases(i,1), test_cases(i,2), test_cases(i,3), max_sig, status);
end


%% 4. WIZUALIZACJA TABELI PRAWDY I STOSUNKU ON/OFF
% Wykres podsumowujący wyniki z Sekcji 3.
figure('Position', [100 100 1400 600], 'Name', 'Test wszystkich kombinacji');

% Przygotowanie danych do wykresu
signals = [results.max_signal];
labels = cell(1, length(results));
colors = zeros(length(results), 3); 

% Pętla do kolorowania słupków 
for i = 1:length(results)
    labels{i} = sprintf('%d%d%d', results(i).miR21, results(i).miR155, results(i).miR34);
    
    % Logika kolorowania poprzez odczytywanie stanu z wyników
    switch results(i).expected_state
        case 1  % Stan ON
            colors(i,:) = [0 0.8 0];  % Zielony
        case 2  % Stan INHIBITED
            colors(i,:) = [1 0 0];    % Czerwony
        case 0  % Stan OFF
            colors(i,:) = [0.3 0.3 0.3];  % Szary
    end
end

% Wykres 1: Sygnał dla wszystkich 8 kombinacji
subplot(1,2,1);
b = bar(1:8, signals);
set(gca, 'YScale', 'log', 'XTick', 1:8, 'XTickLabel', labels);
xlabel('Kombinacja wejść (miR21-miR155-miR34)');
ylabel('Max Signal [M]');
title('Sygnał wyjściowy dla wszystkich kombinacji');
grid on;
ylim([1e-15 1e-8]);

% Aplikowanie zdefiniowanych kolorów do słupków
h = findobj(gca,'Type','bar');
h.FaceColor = 'flat';
h.CData = colors;

% Wykres 2: Obliczenie i wizualizacja stosunku ON/OFF
subplot(1,2,2);

% Obliczenie ratio przy użyciu filtrowania logicznego
all_signals = [results.max_signal];
all_states = [results.expected_state];

% Znalezienie sygnału ON (gdzie stan == 1)
on_signal = mean(all_signals(all_states == 1));

% Znalezienie sygnałów OFF (gdzie stan == 0 lub == 2)
off_signals = all_signals(all_states == 0 | all_states == 2);
avg_off = mean(off_signals);

% Wygenerowanie wykresu słupkowego ON vs OFF
bar([1 2], [on_signal, avg_off]);
set(gca, 'YScale', 'log', 'XTick', [1 2], 'XTickLabel', {'ON (1,1,0)', 'Avg OFF'});
ylabel('Signal [M]');
title(sprintf('ON/OFF Ratio: %.1fx', on_signal/avg_off));
grid on;


%% 5. ANALIZA DYNAMIKI SYSTEMU (STAN ON)
% Bada, jak szybko system reaguje (dynamika czasowa) w pożądanym 
% stanie 'ON' (1,1,0).
fprintf('\nAnaliza czasów odpowiedzi (dla stanu ON)\n');

stop_time_dynamic = 20; 

% Przygotowanie modelu dla stanu ON 
m_on = copyobj(model);
s = m_on.Species;
set(s(strcmp({s.Name}, 'miR21')), 'InitialAmount', CONFIG.conc_high); 
set(s(strcmp({s.Name}, 'miR155')), 'InitialAmount', CONFIG.conc_high);
set(s(strcmp({s.Name}, 'miR34')), 'InitialAmount', CONFIG.conc_low);

% Konfiguracja solwera
cs_on = getconfigset(m_on, 'active');
cs_on.SolverType = 'ode15s';
cs_on.StopTime = stop_time_dynamic; 
cs_on.SolverOptions.AbsoluteTolerance = 1e-20;
cs_on.SolverOptions.RelativeTolerance = 1e-12;
cs_on.SolverOptions.MaxStep = 0.5;

% Uruchomienie symulacji
sim_obj = sbiosimulate(m_on);

% Ekstrakcja danych czasowych i sygnałowych
time = sim_obj.Time;
signal_idx = find(strcmp(sim_obj.DataNames, 'Signal'), 1);
signal_data = sim_obj.Data(:, signal_idx);
signal_data(signal_data < 0) = 0; % Oczyszczenie danych

% Obliczenie kluczowych metryk: t50 i t90 (czas do 50%/90% sygnału)
max_val = max(signal_data);
t50_idx = find(signal_data >= 0.5*max_val, 1, 'first');
t90_idx = find(signal_data >= 0.9*max_val, 1, 'first');

if ~isempty(t50_idx)
    t50 = time(t50_idx);
    fprintf('Czas do 50%% max sygnału (t50):  %.3f s\n', t50);
end
if ~isempty(t90_idx)
    t90 = time(t90_idx);
    fprintf('Czas do 90%% max sygnału (t90):  %.3f s\n', t90);
end

% Generowanie wykresu odpowiedzi czasowej
figure('Position', [150 150 1000 500], 'Name', 'Odpowiedź czasowa');
plot(time, signal_data, 'b-', 'LineWidth', 2); hold on;
% Zaznaczenie punktów t50 i t90 na wykresie
if ~isempty(t50_idx)
    plot(t50, signal_data(t50_idx), 'ro', 'MarkerSize', 10, 'DisplayName', '50% max');
end
if ~isempty(t90_idx)
    plot(t90, signal_data(t90_idx), 'go', 'MarkerSize', 10, 'DisplayName', '90% max');
end
xlabel('Czas [s]');
ylabel('Signal [M]');
title('Odpowiedź czasowa systemu (stan ON: 1,1,0)');
legend('Location', 'best');
grid on;
set(gca, 'YScale', 'log');


%% 6. ANALIZA WRAŻLIWOŚCI 'AND' (EC50)
% Bada, jakie stężenie wejściowe (miR21 i miR155) jest potrzebne do 
% aktywacji systemu. Oblicza EC50 (stężenie dla 50% max odpowiedzi).
fprintf('\nAnaliza wrażliwości na stężenie wejść AND\n');

% Definicja parametrów 
conc_range_ec50 = logspace(-9, -5, 10);  % Od 1nM do 10μM
stop_time_ec50 = 30; 

sensitivity_results = zeros(length(conc_range_ec50), 1);

% Pętla symulacji
for i = 1:length(conc_range_ec50)
    [~, max_sig] = simulate_scenario(model, ...
        conc_range_ec50(i), ...
        conc_range_ec50(i), ...
        CONFIG.conc_low, ... 
        stop_time_ec50);
    
    sensitivity_results(i) = max_sig;
end

% Generowanie wykresu krzywej odpowiedzi
figure('Position', [200 200 800 600], 'Name', 'Wrażliwość (EC50)');
semilogx(conc_range_ec50*1e9, sensitivity_results, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Stężenie wejściowe miR21/miR155 [nM]');
ylabel('Max Signal [M]');
title('Krzywa wrażliwości biokomputera (EC50)');
grid on;
set(gca, 'YScale', 'log');

% Obliczenie EC50
max_response = max(sensitivity_results);
[~, ec50_idx] = min(abs(sensitivity_results - 0.5*max_response));
ec50 = conc_range_ec50(ec50_idx) * 1e9;  % Konwersja na nM
fprintf('Obliczone EC50 (50%% max odpowiedzi): %.2f nM\n', ec50);


%% 7. EKSPORT GŁÓWNYCH WYNIKÓW

% Tworzenie tabeli wyników
T = table();
T.miR21 = [results.miR21]';
T.miR155 = [results.miR155]';
T.miR34 = [results.miR34]';
T.Max_Signal_M = [results.max_signal]';

% Konwersja kodów stanu (0,1,2) na etykiety tekstowe
state_numbers = [results.expected_state]';
state_labels = cell(length(state_numbers), 1);
state_labels(state_numbers == 0) = {'OFF'};
state_labels(state_numbers == 1) = {'ON'};
state_labels(state_numbers == 2) = {'INHIBITED'};
T.Expected_State = state_labels;

% Zapis tabeli do pliku CSV
writetable(T, 'biocomputer_results_summary.csv');

% Zapisywanie wykresu poprzez wyszukiwanie po nazwie)
fig_to_save = findobj('Type', 'figure', 'Name', 'Wrażliwość (EC50)');
if ~isempty(fig_to_save)
    exportgraphics(fig_to_save, 'biocomputer_sensitivity.png', 'Resolution', 300);
else
    warning('Nie znaleziono figury o nazwie "Wrażliwość (EC50)", plik nie został zapisany.');
end

% Użycie zmiennych 'on_signal' i 'avg_off' z Sekcji 4
inhibition_signal = mean([results([results.expected_state] == 2).max_signal]);
fprintf('Finalne ON/OFF ratio: %.1fx\n', on_signal / avg_off);
fprintf('Skuteczność inhibicji przez miR34: %.1f%%\n', (1 - inhibition_signal / on_signal)*100);


%% 8. ANALIZA WRAŻLIWOŚCI 'NOT' (IC50)
% Sprawdzenie, jakie stężenie inhibitora (miR34) jest potrzebne
% do 50% wyłączenia systemu, gdy jest on w stanie ON (1,1,0).

% Definicja parametrów 
conc_range_inhibitor = logspace(-10, -5, 20); % od 100pM do 10μM
stop_time_ic50 = 30; % Czas symulacji
inhibitor_results = zeros(length(conc_range_inhibitor), 1);

% Pętla symulacji inhibicji
for i = 1:length(conc_range_inhibitor)
    
    % Użycie globalnej zmiennej
    [~, max_sig] = simulate_scenario(model, ...
        CONFIG.conc_high, ...
        CONFIG.conc_high, ...
        conc_range_inhibitor(i), ...
        stop_time_ic50);
    
    inhibitor_results(i) = max_sig;
end

% Generowanie wykresu krzywej inhibicji
figure('Position', [250 250 800 600], 'Name', 'Wrażliwość na Inhibitor (IC50)');
semilogx(conc_range_inhibitor * 1e9, inhibitor_results, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Stężenie inhibitora [miR34] (nM)');
ylabel('Max Signal [M]');
title('Krzywa Inhibicji Biokomputera (IC50)');
grid on;
set(gca, 'YScale', 'log');

% Obliczenie IC50
max_response = max(inhibitor_results);
[~, ic50_idx] = min(abs(inhibitor_results - 0.5*max_response));
ic50 = conc_range_inhibitor(ic50_idx) * 1e9; % Konwersja na nM
fprintf('Obliczone IC50 (50%% inhibicji): %.2f nM\n', ic50);

% Zapis wykresu
exportgraphics(gcf, 'biocomputer_inhibition_IC50.png', 'Resolution', 300);

%% 9. ANALIZA ROBUSTNOŚCI (STRESS TEST) - WPŁYW PRZECIEKU
% Badanie, jak parametr 'k_leak_gate21' (błędne generowanie sygnału)
% wpływa na jakość systemu (stosunek ON/OFF).

% Definicja zakresu testu dla parametru 'k_leak_gate21'
leak_range = logspace(-8, -2, 20);
results_on_signal = zeros(length(leak_range), 1);
results_avg_off = zeros(length(leak_range), 1);
results_ratio = zeros(length(leak_range), 1);

% Pętla "Stress Testu"
for i = 1:length(leak_range)
    current_leak_val = leak_range(i);
    
    % Kopia modelu do modyfikacji
    m_leaky = copyobj(model); 
    
    % Modyfikacja parametru 'k_leak_gate21' w kopii modelu
    param_obj = sbioselect(m_leaky, 'Type', 'parameter', 'Name', 'k_leak_gate21');
    if isempty(param_obj)
        error('Nie można znaleźć parametru "k_leak_gate21" w modelu.');
    end
    set(param_obj, 'Value', current_leak_val);
    
    % Symulacja stanu ON 
    [~, on_signal] = simulate_scenario(m_leaky, ...
        CONFIG.conc_high, ...
        CONFIG.conc_high, ...
        CONFIG.conc_low, ...
        CONFIG.stop_time);
    
    % Definicja 7 stanów OFF
    test_cases_off = [
        0 0 0; 0 0 1; 0 1 0; 0 1 1;
        1 0 0; 1 0 1; 1 1 1
    ];
    off_signals = zeros(7, 1);
    % Pętla symulująca wszystkie 7 stanów OFF
    for j = 1:7
        miR21_val = test_cases_off(j,1) * CONFIG.conc_high;
        miR155_val = test_cases_off(j,2) * CONFIG.conc_high;
        miR34_val = test_cases_off(j,3) * CONFIG.conc_high;
        
        [~, max_sig] = simulate_scenario(m_leaky, miR21_val, miR155_val, miR34_val, CONFIG.stop_time);
        off_signals(j) = max_sig;
    end
    
    avg_off = mean(off_signals);
    
    % Zapis wyników pętli
    results_on_signal(i) = on_signal;
    results_avg_off(i) = avg_off;
    results_ratio(i) = on_signal / avg_off;
end

% Generowanie podwójnego wykresu
figure('Position', [300 300 1000 600], 'Name', 'Analiza robustności - Przeciek (Leakage)');

% Wykres: Sygnał ON vs Sygnał OFF
subplot(1, 2, 1);
loglog(leak_range, results_on_signal, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(leak_range, results_avg_off, 'r-x', 'LineWidth', 2, 'MarkerSize', 8);
hold off;
xlabel('Wartość parametru przecieku [k\_leak\_gate21]');
ylabel('Max Signal [M]');
title('Wpływ przecieku na sygnał ON vs OFF');
legend('Sygnał ON (1,1,0)', 'Średni sygnał OFF', 'Location', 'northwest');
grid on;

% Wykres: Stosunek ON/OFF
subplot(1, 2, 2);
loglog(leak_range, results_ratio, 'b-d', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Wartość parametru przecieku [k\_leak\_gate21]');
ylabel('Stosunek ON/OFF Ratio');
title('Wpływ przecieku na ON/OFF Ratio');
grid on;

% Zapis wykresu
exportgraphics(gcf, 'biocomputer_robustness_LEAK.png', 'Resolution', 300);


%% 10. ANALIZA ROBUSTNOŚCI - WPŁYW CROSSTALK
% Badanie, jak kradzież wejścia miR21 przez reakcje 
% uboczne ('k_cross') wpływa na zdolność systemu do aktywacji.

% Definicja zakresu testu dla parametru 'k_cross'
cross_range = logspace(1, 8, 20); 
results_on_signal_c = zeros(length(cross_range), 1);
results_avg_off_c = zeros(length(cross_range), 1);
results_ratio_c = zeros(length(cross_range), 1);

% Pętla "Stress Testu"
for i = 1:length(cross_range)
    current_cross_val = cross_range(i);
    
    % Kopia modelu do modyfikacji
    m_cross = copyobj(model); 
    
    % Modyfikacja parametru 'k_cross' w kopii modelu
    param_obj = sbioselect(m_cross, 'Type', 'parameter', 'Name', 'k_cross_21_155');
    if isempty(param_obj)
        error('Nie można znaleźć parametru "k_cross" w modelu.');
    end
    set(param_obj, 'Value', current_cross_val);
    
    % Symulacja stanu ON 
    [~, on_signal] = simulate_scenario(m_cross, ...
        CONFIG.conc_high, ...
        CONFIG.conc_high, ...
        CONFIG.conc_low, ...
        CONFIG.stop_time);
    
    % Definicja 7 stanów OFF
    test_cases_off = [
        0 0 0; 0 0 1; 0 1 0; 0 1 1;
        1 0 0; 1 0 1; 1 1 1
    ];
    off_signals = zeros(7, 1);
    % Pętla symulująca wszystkie 7 stanów OFF
    for j = 1:7
        miR21_val = test_cases_off(j,1) * CONFIG.conc_high;
        miR155_val = test_cases_off(j,2) * CONFIG.conc_high;
        miR34_val = test_cases_off(j,3) * CONFIG.conc_high;
        
        [~, max_sig] = simulate_scenario(m_cross, miR21_val, miR155_val, miR34_val, CONFIG.stop_time);
        off_signals(j) = max_sig;
    end
    
    avg_off = mean(off_signals);
    
    % Zapis wyników pętli
    results_on_signal_c(i) = on_signal;
    results_avg_off_c(i) = avg_off;
    results_ratio_c(i) = on_signal / avg_off;
end

% Generowanie podwójnego wykresu
figure('Position', [350 350 1000 600], 'Name', 'Analiza robustności - Crosstalk');

% Wykres: Sygnał ON vs Sygnał OFF
subplot(1, 2, 1);
loglog(cross_range, results_on_signal_c, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
loglog(cross_range, results_avg_off_c, 'r-x', 'LineWidth', 2, 'MarkerSize', 8);
hold off;
xlabel('Wartość parametru crosstalk [k\_cross]');
ylabel('Max Signal [M]');
title('Wpływ crosstalk na sygnał ON vs OFF');
legend('Sygnał ON (1,1,0)', 'Średni sygnał OFF', 'Location', 'southwest');
grid on;

% Wykres: Stosunek ON/OFF
subplot(1, 2, 2);
loglog(cross_range, results_ratio_c, 'b-d', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Wartość parametru crosstalk [k\_cross]');
ylabel('Stosunek ON/OFF Ratio');
title('Wpływ crosstalk na ON/OFF Ratio');
grid on;

% Zapis wykresu
exportgraphics(gcf, 'biocomputer_robustness_CROSSTALK.png', 'Resolution', 300);