//
// Created by Benedikt Völker on 25.04.17.
//

#include <queue>
#include "CMyMath.h"
#include "Risiko_Spiel.h"
#define RED     "\033[31m"
#define RESET   "\033[0m"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

CMyMath::CMyMath() {
    time_t t;
    time(&t);
    srand((unsigned int)t);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

string CMyMath::cs ( double zahl ){
    string res = to_string(zahl);

    for(char& c : res) {
        if ( c == '.')
            c = ',';
    }

    return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::integral_e ( double x ){
    double res = 1;

    for(int i = 1; i < 20; i++){
        res += pow(x, i) / double( fakultaet(i) );
    }

    return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::F ( double x ){

    return integral_e(-(x*x/2));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::NueberK( int n, int k ){

    double zaehler = 1;
    double nenner = 1;

    for(;k > 0; k--, n--){
        zaehler *= n;
        nenner *= k;
    }

    return zaehler / nenner;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CMyMath::fakultaet(int n) {
    if (n == 0)
        return 1;
    return n * fakultaet(n - 1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CMyMath::natuerlicheZahlinN ( int n ) {
    return rand() % n + 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CMyMath::natuerlicheZahlinN_0 ( int n ) {
    return rand() % (n+1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CMyMath::natuerlicheZahlinN_neu ( int n ) {
    int range = RAND_MAX / n;
    int i = 1;
    int random = rand();

    while( i * range < random ){
        i++;
        if(i > n){
            i = 1;
            random = rand();
        }
    }
    return i;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CMyMath::natuerlicheZahlinAB ( int a, int b ) {
    int res = 0;
    while( !( res >= a && res <= b ) ){
        res = natuerlicheZahlinN( b );
    }
    return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::zufallszahl_in_0_1 (){
    return double(rand()) / double(RAND_MAX);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool CMyMath::bernoulliSimulation ( double p ) {
    return static_cast <double> (rand()) / static_cast <double> (RAND_MAX) <= p;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::generiereExponentialverteilung ( double lambda ){
    return - log( zufallszahl_in_0_1() ) / lambda;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::generiereNormalverteilung ( double mittelwert, double varianz ){
    return mittelwert + sqrt( -2 * varianz * log( zufallszahl_in_0_1() ) ) * sin( 2 * M_PI * zufallszahl_in_0_1() );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::generiereHistogrammExponentialverteilt ( int iterationen, int a, int b, int anz, double lambda ){
    double step = (b-a)/double(anz);
    double table[anz][2];

    for(int i = 0; i < anz; i++){
        table[i][0] = a + (i+1) * step;
        table[i][1] = 0;
    }

    double zahl;
    for(int i = 0; i < iterationen; i++){
        zahl = generiereExponentialverteilung( lambda );
        for(int j = 0; j < anz; j++){
            if(table[j][0] > zahl ){
                table[j][1] += 1;
                break;
            }
        }
    }

    double scale = 1/table[0][1] * lambda;

    for(int i = 0; i < anz; i++){
        // cout << "Zwischen " << table[i][0] - step << " und " << table[i][0] << " : " << (table[i][1] * scale) << endl;
        cout << cs(table[i][1] * scale) << endl;
    }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::generiereHistogrammNormalverteilt ( int iterationen, int dx, int anz, int mittel, double varianz){
    double step = (2*dx)/double(anz);
    double table[50][2];

    for(int i = 0; i < anz; i++){
        table[i][0] = - dx + (i+1) * step;
        table[i][1] = 0;
    }

    double zahl;
    for(int i = 0; i < iterationen; i++){
        zahl = generiereNormalverteilung(mittel,varianz);
        for(int j = 0; j < anz; j++){
            if(zahl > - dx && table[j][0] > zahl ){
                table[j][1] += 1;
                break;
            }
        }
    }

    int point = (int)(50 * ( dx+mittel ) / ( 2*dx ));

    double scale = 1/table[(int)(point)][1] / (sqrt(2*M_PI*varianz));

    for(int i = 0; i < anz; i++){
        // cout << "Zwischen " << table[i][0] - step << " und " << table[i][0] << " : " << (table[i][1] * scale) << endl;
        cout << cs(table[i][1] * scale) << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::simuliereOptimaleWartezeitExponentiell (int b, double lambda, int bereich, double genauigkeit, int iterationen ){
    int anz = bereich / genauigkeit;
    double table[anz][2];

    for(int i = 0; i < anz; i++){
        table[i][0] = - bereich + (i+1) * genauigkeit;
        table[i][1] = 0;
    }

    double x;
    for(int i = 0; i < anz; i++){
        x = table[i][0];
        for(int j = 0; j < iterationen; j++){
            table[i][1] += wartezeitExponentiell(b, x, lambda);
        }
        table [i][1] /= iterationen;
    }

    int minimum = 0;
    for(int i = 0; i < anz; i++){
        // cout << "Zwischen " << table[i][0] - step << " und " << table[i][0] << " : " << (table[i][1] * scale) << endl;
        if(table[i][1] < table[minimum][1]){
            minimum = i;
        }
        cout << cs(table[i][1]) << endl;
    }

    cout << "Die minimale Wartezeit liegt bei " << - table[minimum][0] << " Minuten vor Abfahrt." << endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::wartezeitExponentiell(int b, double x, double lambda){
    double deltaAnkunft = fmod(x + generiereExponentialverteilung( lambda ), b);
    if (deltaAnkunft <= 0){
        return - deltaAnkunft;
    }
    else {
        return b - deltaAnkunft;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::simuliereOptimaleWartezeitNormalverteilt(int b, double stdAbweichung, int bereich, double genauigkeit, int iterationen) {
    int anz = bereich / genauigkeit;
    double table[anz][2];

    for(int i = 0; i < anz; i++){
        table[i][0] = - bereich + (i+1) * genauigkeit;
        table[i][1] = 0;
    }

    double x;
    for(int i = 0; i < anz; i++){
        x = table[i][0];
        for(int j = 0; j < iterationen; j++){
            table[i][1] += wartezeitNormalverteilt(b, x, stdAbweichung);
        }
        table [i][1] /= iterationen;
    }

    int minimum = 0;
    for(int i = 0; i < anz; i++){
        // cout << "Zwischen " << table[i][0] - step << " und " << table[i][0] << " : " << (table[i][1] * scale) << endl;
        if(table[i][1] < table[minimum][1]){
            minimum = i;
        }
        cout << cs(table[i][1]) << endl;
    }

    cout << "Die minimale Wartezeit liegt bei " << cs(- table[minimum][0]) << " Minuten vor Abfahrt und beträgt " << cs(table[minimum][1]) << " Minuten." << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::wartezeitNormalverteilt(int b, double x, double stdAbweichung) {
    double deltaAnkunft = generiereNormalverteilung( x, stdAbweichung*stdAbweichung );
    if (deltaAnkunft <= 0){
        return - deltaAnkunft;
    }
    else {
        return b - deltaAnkunft;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// PRAKTIKUM 1 : Aufgabe 1

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::wahrscheinlichkeitRisiko(int k, int l) {

    int teilsiege, moeglichkeiten = 0;

    vector<int> angreifer = vector<int>(k,1);
    vector<int> verteidiger = vector<int>(l,1);
    vector<int> results = vector<int>(l+1);
    vector<int> ang_sort;
    vector<int> ver_sort;

    while(angreifer.at(k-1) < 7) {
        bool naechster_angriff = false;
        while (!naechster_angriff) {
            moeglichkeiten++;

            ang_sort = angreifer;
            sort(ang_sort.begin(), ang_sort.end(), std::greater<int>());

            ver_sort = verteidiger;
            sort(ver_sort.begin(), ver_sort.end(), std::greater<int>());

            teilsiege = 0;

            for(int i = 0; i < l; i++){
                if(ang_sort.at(i) > ver_sort.at(i)){
                    teilsiege++;
                }
            }

            results.at(teilsiege)++;

            verteidiger.at(0)++;

            for(int i = 0; i < l; i++){
                if(verteidiger.at(i) > 6){
                    verteidiger.at(i) = 1;
                    if(i+1 < l) {
                        verteidiger.at(i+1)++;
                    }
                    else {
                        naechster_angriff = true;
                    }
                }
            }
        }

        angreifer.at(0)++;

        for(int i = 0; i < k; i++){
            if( i+1 < k && angreifer.at(i) > 6 ){
                angreifer.at(i) = 1;
                angreifer.at(i+1)++;
            }
        }
    }

    for(int i = 0; i <= l; i++){
        cout << "Wahrscheinlichkeit für einen "<< i << "-fach Sieg: " << results.at(i)/double(moeglichkeiten) << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::wahrscheinlichkeitBitfehlerN32( double p, int f, int n ){
    int unbrauchbar = 0;
    int bitfehler = 0;

    for(int i = 0; i < n; i++){
        for(int b = 0; b < 32; b++){
            if(bernoulliSimulation(p))
                bitfehler++;
        }
        if(bitfehler > f){
            unbrauchbar++;
        }
        bitfehler = 0;
    }

    // cout << "Die Wahrscheinlichkeit für mehr als 2 Bitfehler ist: " << unbrauchbar / double(n) * 100 <<  "%" << endl;
    return unbrauchbar / double(n);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::streuungBitfehler( double p, int f, int n, int N, bool ausgabe){
    vector<double> stichproben = vector<double>(N, 0);
    double mittelwert = 0, varianz = 0, test = 0;

    for(int i = 0; i < N; i++){
        test = wahrscheinlichkeitBitfehlerN32( p, f, n );
        stichproben[i] = test;
        mittelwert += test;
        varianz += test * test;
    }

    mittelwert /= N;

    varianz = 1.0/(N-1.0) * (varianz - N * mittelwert * mittelwert);

    if(ausgabe){
        cout << "Mittelwert: " << mittelwert << endl;
        cout << "Varianz: " << varianz << endl;
        cout << "Die Standartabweichung bei " << n << " durchläufen beträgt: " << sqrt(varianz) << "" << endl;
    }
    else{
        cout << cs(sqrt(varianz)/n) << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::streuungRisiko( int k, int l, int n, int N, bool ausgabe){
    vector<double> stichproben = vector<double>(N, 0);
    double mittelwert = 0, varianz = 0, test = 0;
    Risiko_Spiel spiel = Risiko_Spiel();

    for(int i = 0; i < N; i++){
        test = spiel.spielen(k, l, n, false);
        stichproben[i] = test;
        mittelwert += test;
        varianz += test * test;
    }


    mittelwert /= N;

    varianz = 1.0/(N-1.0) * (varianz - N * mittelwert * mittelwert);

    if(ausgabe){
        cout << "Mittelwert: " << mittelwert << endl;
        cout << "Varianz: " << varianz << endl;
        cout << "Die Standartabweichung bei " << n << " durchläufen beträgt: " << sqrt(varianz) << "" << endl;
    }
    else{
        cout << cs(sqrt(varianz)/n) << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

mittelwert_varianz CMyMath::stichprobeNormalverteilung(int n, double erwartungswert, double standardabweichug ){
    mittelwert_varianz result;
    double beobachtung;

    for(int i = 0; i < n; i++){
        beobachtung = generiereNormalverteilung(erwartungswert, standardabweichug * standardabweichug);
        result.mittelwert += beobachtung;
        result.varianz += beobachtung * beobachtung;
    }

    result.mittelwert /= n;
    result.varianz = 1.0/(n-1.0) * ( result.varianz - n * result.mittelwert * result.mittelwert);

    return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::stichprobemMittelNormalverteilung(int N, int n, double erwartungswert, double standardabweichung ){
    double mittelwert = 0, std_abw = 0;
    mittelwert_varianz result;

    for(int i = 0; i < N; i++){
        result = stichprobeNormalverteilung( n, erwartungswert, standardabweichung );
        mittelwert += result.mittelwert;
        std_abw += sqrt( result.varianz );
    }

    mittelwert /= N;
    std_abw /= N;

    cout << "Mittelwert:" << endl
         << "Original: " << erwartungswert << endl
         << "Erreicht: " << mittelwert << endl
         << "Differenz: " << erwartungswert - mittelwert << endl << endl
         << "Varianz:" << endl
         << "Original: " << standardabweichung << endl
         << "Erreicht: " << std_abw << endl
         << "Differenz: " << standardabweichung - std_abw << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::stichprobeExponentialverteilung ( int n, double lambda ){
    double mittelwert = 0, c = (double(n)-1)/double(n);

    for(int i = 0; i < n; i++){
        mittelwert += generiereExponentialverteilung( lambda );
    }

    mittelwert /= n;

    return c/mittelwert;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::stichprobemMittelExponentialverteilung (int N, int n, double lambda ){
    double mittelwert = 0;

    for(int i = 0; i < N; i++){
        mittelwert += stichprobeExponentialverteilung( n, lambda );
    }

    mittelwert /= N;

    cout << "Mittelwert:" << endl
         << "Original: " << lambda << endl
         << "Erreicht: " << mittelwert << endl
         << "Differenz: " << lambda - mittelwert << endl;
    cout << "Prozentuale Abweichung: " << (lambda - mittelwert)/mittelwert << endl;
    cout << "lambda / mittelwert : " << lambda/mittelwert << endl;
    cout << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::simuliereKonfidenzbereich( int N, int n ) {
    int K = 0;
    double p = 0;

    for(int i = 0; i < N; i++ ) {
        p = zufallszahl_in_0_1();
        if( stichprobeKonfidenzbereich( n, p ) ) {
            K++;
        }
    }

    double p_pinKonf = double(K)/double(N) * 100;

    if( p_pinKonf < 95 ) {
        cout << RED << "Erfolgswahrscheinlichkeit 'p' lag zu " << p_pinKonf << "% im Konfidenzbereich." << RESET << endl;
    }
    else {
        cout << "Erfolgswahrscheinlichkeit 'p' lag zu " << p_pinKonf << "% im Konfidenzbereich." << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool CMyMath::stichprobeKonfidenzbereich( int n, double p ) {
    int k = 0;

    for (int i = 0; i < n; i++) {
        if (bernoulliSimulation(p)) {
            k++;
        }
    }

    konfidenzbereich p12 = berechneKonfidenzbereich( n, k );

    //cout << "Konfidenzbereich von p( " << p << " ): [ " << p12.begin << "; " << p12.end << " ]" << endl;

    return p12.begin < p && p < p12.end;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

konfidenzbereich CMyMath::berechneKonfidenzbereich( int n, int k ) {
    double c = 1.96; // Qantil( 1 - a/2 = 0.975 )

    konfidenzbereich p12;
    if( k < 50 || n-k < 50 ) {
        double x = 1.0 / (c * c + double(n));
        double y = double(k) + (c * c / 2.0);
        double z = c * sqrt(double(k) * double(n - k) / double(n) + c * c / 4.0);

        p12.begin = x * (y - z);
        p12.end = x * (y + z);
    } else {
        double x = 0, y = 0, z = 0;
        if( true ) {
            x = double(k) / double(n);
            y = c / sqrt(n);
            z = sqrt(x * (1 - x));
        } else {
            x = 1.0 / double(n);
            y = k;
            z = c * sqrt(double(k) * double(n - k) / double(n));
        }
        p12.begin = x - y*z;
        p12.end = x + y*z;
    }

    return p12;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

akzeptanzbereich CMyMath::berechneAkzeptanzbereich( int n, double p, double a, int mode ) {

    akzeptanzbereich AK = akzeptanzbereich();
    //double k0_calc = 0;
    //double k1_calc = 0;

    switch( mode ){
        case -1:
            AK.k0 = berechneK0( n, p, 1-a );
            AK.k1 = n;
            break;
        case 1:
            AK.k0 = 0;
            AK.k1 = berechneK1( n, p, 1-a);
            break;
        default:
            AK.k0 = berechneK0( n, p, (1-a)/2 );
            AK.k1 = berechneK1( n, p, (1-a)/2 );
            break;
    }

    // Berechnung Akzeptanzbereich

    double k0_calc = double(n) * p - ( sqrt( double(n) * p * (1-p) ) * 1.645 + 0.5);
    double k1_calc = double(n) * p + ( sqrt( double(n) * p * (1-p) ) * 1.645 + 0.5);
    /*
    cout << "K0: " << AK.k0 << endl
         << "K1: " << AK.k1 << endl
         << "K0_calc: " << k0_calc << endl
         << "K1_calc: " << k1_calc << endl;
    */
    return AK;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::berechneK0( int n, double p, double a ) {
    int i = 0;
    double k0 = NueberK(n, i) * pow(p, i) * pow(1-p, n-i);

    while(k0 < a){
        i++;
        k0 += NueberK(n, i) * pow(p, i) * pow(1-p, n-i);
    }

    return double(i);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double CMyMath::berechneK1( int n, double p, double a ) {
    int i = n;
    double k1 = 1 - NueberK(n, i) * pow(p, i) * pow(1-p, n-i);

    while(k1 > 1-a){
        i--;
        k1 -= NueberK(n, i) * pow(p, i) * pow(1-p, n-i);
    }

    return double(i);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool CMyMath::hypothesentest() {
    double p = zufallszahl_in_0_1();
    double a = 0.9; //natuerlicheZahlinAB( 10, 100);
    //a = 1 - 1/a;
    int n = natuerlicheZahlinAB( 20, 200 );
    akzeptanzbereich AB = berechneAkzeptanzbereich( n, p, a, 0 );

    /*
    cout << endl << "Generierte Daten: " << endl
         << "p: " << p << endl
         << "a: " << a << endl
         << "n: " << n << endl
         << "k0: " << AB.k0 << endl
         << "k2: " << AB.k1 << endl
         << endl;
    */

    int k = 0;

    for(int i = 0; i < n; i++){
        k += int(bernoulliSimulation( p ));
    }

    double p_sim = double(k) / double(n);
    /*
    cout << "k: " << k << endl
         << "n: " << n << endl
         << "Simuliertes p: " << p_sim << endl
         << endl;
    */

    if( AB.k0 <= k && k <= AB.k1 ){
        //cout << "K liegt im Konfidenzbereich." << endl << endl;
        return true;
    } else {
        //cout << "K liegt NICHT im Konfidenzbereich." << endl << endl;
        return false;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int CMyMath::generierePoisson(double a, double t) {
    double p = zufallszahl_in_0_1();
    double d = a * t;
    int k = 0;
    double poisson = pow(d,k) / double( fakultaet( k ) )  * exp(-d);;

    while(p > poisson ){
        k++;
        poisson += pow(d,k) / double( fakultaet( k ) )  * exp(-d);
    }

    return k;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CMyMath::mittlereAufenthaltszeit( int N, double a, double b, int mode) {

    queue<kunde> fifo = queue<kunde>();

    int n_come = 0;
    int n_serv = 0;

    int counter = 0;
    int c2 = 0;

    double d_wartezeit = 0;

    for(double t = 0; t < N; t = t + 0.1){
        n_come = generierePoisson( a , 0.1);

        switch( mode ){
            case 0:
                n_serv = generierePoisson( b , 0.1);
                break;
            case 1:
                if( t - int(t) < 0.01 )
                    n_serv = natuerlicheZahlinN_0( 2 / b );
                else
                    n_serv = 0;
                break;
                break;
            case 2:
                //cout << t << " - " << (int)t << " = " << t - int(t) << endl;
                if( c2 == 10 ) {
                    //cout << t << endl;
                    c2 = 0;
                    n_serv = b;
                } else {
                    n_serv = 0;
                    c2++;
                    //cout << "." ;
                }
                break;
        }


        for(int j = 0; j < n_come; j++){
            kunde k_neu;
            k_neu.t0 = t;
            fifo.push( k_neu );
        }

        for(int j = 0; j < n_serv && fifo.size() > 0 ; j++){
            //cout << fifo.front().t0 << " - " << t << " = " << t - fifo.front().t0 << endl;
            if( t - fifo.front().t0 > 1000)
                cout << endl;
            d_wartezeit += t - fifo.front().t0;
            fifo.pop();
            counter++;
        }

        //cout << fifo.size() << endl;
    }

    d_wartezeit /= double(counter);

    cout << "Durchschnittliche Wartezeit: " << d_wartezeit << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////