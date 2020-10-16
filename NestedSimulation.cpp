


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

using namespace std;

class BlackScholes
{
public:
    BlackScholes(double SigneInput, double SInput, double KInput, double rInput, double volInput, double TInput);
    BlackScholes();
    ~BlackScholes();
    double d1();
    double d2();
    double Prix();
    void PrTest();
    void Erreur();
    double Signe;
private:
    double S, K, r, vol, T;
};




BlackScholes::BlackScholes(double SigneInput, double SInput, double KInput, double rInput, double volInput, double TInput)//constructeur
{
    Signe = SigneInput;
    S = SInput;
    K = KInput;
    r = rInput;
    vol = volInput;
    T = TInput;
}

//BlackScholes::BlackScholes() = { 0,0,0,0,0 }//Constructeur par défaut
//{
//}

double BlackScholes::d1()//Calcul de  d1
{
    return (log(S / K) + (r + pow(vol, 2) / 2) * T) / (vol * sqrt(T));
}

double BlackScholes::d2()//calcul de d2
{

    return d1() - vol * sqrt(T);
}

double BlackScholes::Prix()//pricing de l'option
{

    double x = S * N(Signe * d1());
    double y = K * exp(-1 * r * T) * N(Signe * d2());
    if (Signe == 1)
    {
        return Signe * (x - y);
    }
    if (Signe == -1)
    {
        return Signe * (x - y);
    }
}

void BlackScholes::Erreur()
{
    if (S <= 0)
    {
        cout << "Erreur !" << endl;
        exit(1);
    }
    if (K <= 0)
    {
        cout << "Erreur!" << endl;
        exit(1);
    }
    if (vol <= 0)
    {
        cout << "Erreur!" << endl;
        exit(1);
    }
    if (T <= 0)
    {
        cout << "Erreur! " << endl;
        exit(1);
    }
    if (Signe != -1 && Signe != 1)
    {
        cout << "Erreur!" << endl;
        exit(1);
    }
}

void BlackScholes::PrTest()//Permet de vérifier les variables entrées
{
    if (Signe == 1)
    {
        cout << "L'option est un  call." << endl;
    }
    if (Signe == -1)
    {
        cout << "L'option est un  put." << endl;
    }
    cout << "La volatilite est " << vol << endl;
    cout << "Le prix du Spot est " << S << endl;
    cout << "Le prix du Strike est " << K << endl;
    cout << "LIBOR est " << r << endl;
    cout << "La maturite en fonction du nombre d'annee " << T << endl;
    cout << "d1 est egal à  " << d1() << endl;
    cout << "d2 est egal à " << d2() << endl;
    cout << "Le d1 normalise est egal a " << N(Signe * d1()) << endl;
    cout << "Le d2 normalise est egal a " << N(Signe * d2()) << endl;

}

BlackScholes::~BlackScholes()//Destructeur
{

}

double N(double z)//Calcul de la fonction de répartition
{
    const double b1 = 0.319381530;
    const double b2 = -0.356563782;
    const double b3 = 1.781477937;
    const double b4 = -1.821255978;
    const double b5 = 1.330274429;
    const double p = 0.2316419;
    const double c2 = 0.39894228;

    double a = fabs(z);
    double t = 1.0 / (1.0 + a * p);
    double b = c2 * exp((-z) * (z / 2.0));
    double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
    n = 1 - b * n;
    if (z < 0) n = 1 - n;
    return n;
}



double LoiUniform()
{
    return ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
}


double LoiNorm()
{
    double u1 = LoiUniform();
    double u2 = LoiUniform();
    return cos(8. * atan(1.) * u2) * sqrt(-2. * log(u1));
}


double GenerSpot(double mu, double vol, double t)
{
    double spot = 0;
    return spot = exp((mu - (vol * vol) / 2) * t + vol * sqrt(t) * LoiNorm());
}


int  Argmin(vector<double> T)
{
    int min = 0;
    for (int i = 0; i < T.size(); i++)
    {
        if (T[i] < T[min]) { min = i; }
    }
    min = min + 1;
    return min;
}

//double Sum(vector <double> T) { for (int i = 0; i < T.size(); i++) { T[i + 1] += T[i]; } return T[T.size() - 1]; }
//int Sum(vector<int> T) {  for (int i = 0; i < T.size(); i++) { T[i + 1] += T[i]; } return T[T.size() - 1]; }


//double EstimUnif(int N, int M, double Prime, double SpotZero, double K, double c, double Signe, double tau, double T, double mu, double r, double vol)
//{
//    double MCGlobalEstimator = 0;
//    double l = T - tau;
//    double SpotTau = 0;
//    double SpotT = 0;
//    vector<double> Loss(N);
//    for (int i = 0; i < N; i++) { Loss[i] = 0; }
//    for (int i = 0; i < N; i++) {
//
//        SpotTau = SpotZero * GenerSpot(mu, vol, tau);
//        for (int j = 0; j < M; j++) {
//
//            SpotT = SpotTau * GenerSpot(r, vol, l);
//
//        }
//        if ((Signe * K < Signe * SpotT)) { Loss[i] += Prime - Signe * (SpotT - K); }
//    }
//
//    for (int i = 0; i < N; i++) {
//        Loss[i] = (1. / (double)M) * Loss[i];//cout<<Loss[i]<<endl;
//        if (Loss[i] <= c)
//        {
//            ++MCGlobalEstimator;
//        }
//    }
//    return MCGlobalEstimator *= (1. / (double)N);
//
//}
//
//
//double EstimSequentiel( int N, int M0, int Mbarre, double Prime, double SpotZero, double K, double c, double Signe, double tau, double T, double mu, double r, double vol)
//{
//    double MCGlobalEstimator = 0;
//    double SpotTau = 0;
//    double l = T - tau;
//    double SpotT = 0;
//    vector<double> Loss(N);
//    vector <int> M(N);
//    vector<double> ErrMarg(N);
//    vector<double> StanDeviation(N);
//    vector<double> SD(N);
//    for (int i = 0; i < N; i++) { Loss[i] = 0; ErrMarg[i] = 0; StanDeviation[i] = 0; SD[i] = 0; }
//
//    for (int i = 0; i < N; i++) {
//
//        SpotTau = SpotZero * GenerSpot(mu, vol, tau);
//
//        for (int j = 0; j < M0; j++) {
//
//            SpotT = SpotTau * GenerSpot(r, vol, l);
//        }
//        M[i] = M0;
//        if ((Signe * K < Signe * SpotT)) {
//            Loss[i] += (Prime - Signe * (SpotT - K));
//            Loss[i] *= (1. / (double)M[i]);//cout<<"LOSS "<<Loss[i]<<endl;
//            SD[i] = pow(Signe * (SpotT - K) - Loss[i], 2);
//        }
//        StanDeviation[i] = (double)(1 / (M[i] - 1)) * pow(Sum(SD, M[i]), 0.5);
//        ErrMarg[i] = (M[i] / StanDeviation[i]) * fabs(Loss[i] - c); //cout<<"ERREUR "<<ErrMarg[i]<<endl;
//    }
//    while (Sum(M, N) < Mbarre * N) {
//
//        int j = Argmin(ErrMarg, N);
//        double z = SpotZero * GenerSpot(0.08, 0.2, 1) * GenerSpot(0.03, 0.2, 1);
//        if ((Signe * z > Signe * K)) { Loss[j] = (1 / (M[j] + 1)) * (Prime - Signe * (z - K)) + (M[j] / (M[j] + 1)) * Loss[j]; }
//        M[j] = M[j] + 1;
//
//
//        for (int i = 0; i < N; i++)
//        {
//            if (Loss[i] <= c)
//            {
//                ++MCGlobalEstimator;
//            }
//        }
//        return MCGlobalEstimator *= (1. / (double)N);
//    }
//}
//double EstimThres( int N, double Y, double Prime, double SpotZero, double K, double c, double Signe, double tau, double T, double mu, double r, double vol)
//{
//    double MCGlobalEstimator = 0;
//    double SpotTau = 0;
//    double SpotT = 0;
//    vector<double> Loss(N);
//    vector<int> M(N);
//    vector<double> ErrMarg(N);
//    vector<double> StanDeviation(N);
//    vector<double> SD(N);
//    double l = T - tau;
//    for (int i = 0; i < N; i++) { M[i] = 0; Loss[i] = 0; ErrMarg[i] = 0; StanDeviation[i] = 0; SD[i] = 0; }
//
//    for (int i = 0; i < N; i++) {
//
//        SpotTau = SpotZero * GenerSpot(mu, vol, tau);
//
//        SD[i] = pow((Signe * (SpotTau - K)) - Loss[i], 2);
//        StanDeviation[i] = (double)(1 / (M[i] - 1)) * pow(Sum(SD, M[i]), 0.5);
//        ErrMarg[i] = (M[i] / StanDeviation[i]) * fabs(Loss[i] - c);
//        M[i] = 0;
//
//        do {
//            SpotT = SpotTau * GenerSpot(r, vol, l);
//            M[i] = M[i] + 1;
//            if ((Signe * K < Signe * SpotT)) { Loss[i] += (Prime - Signe * (SpotT - K)); }
//            Loss[i] *= (1. / (double)M[i]);
//
//        } while (ErrMarg[i] < Y);
//    }
//    for (int i = 0; i < N; i++)
//    {
//        if (Loss[i] <= c)
//        {
//            ++MCGlobalEstimator;
//        }
//    }
//    return MCGlobalEstimator *= (1. / (double)N);
//
//}

int main() {


    double Signe, S, K, r, vol, T;//, c, tau, mu;
    cout << "Suivez les instructions qui vont vous permettre de calculer d'abord" << endl << "la valeur de votre portefeuille a l'instant 0 " << endl << "qui sera constitue soit d'un Put ou d'un Call" << endl << "puis evaluer la proba de perte pour le seuil c" << endl << "de votre portefeuille a l'instant tau " << endl;
    cout << endl << endl;
    cout << "Faites entrer -1 si votre portefeuille contient un Put" << endl << "ou 1 s'il contient un Call" << endl;
    cin >> Signe;
    cout << "Faites entrer le prix du Spot" << endl;
    cin >> S;
    cout << "Faites entrer le prix du Strike" << endl;
    cin >> K;
    cout << "Faites entrer le taux Libor  (Si le Libor est de 4%,  entrez .04)" << endl;
    cin >> r;
    cout << "Faites entrez la Vol (si la vol est de 20% , entrez .2)" << endl;
    cin >> vol;
    cout << "Faites entrez la Maturite de l'option en annees " << endl;
    cin >> T;

    cout << endl;
    BlackScholes Test(Signe, S, K, r, vol, T);
    Test.Erreur();
    double Prime = Test.Prix();
    cout << endl;
    cout << "La valeur initiale de votre portefeuille est egal a " << Prime << " euros" << endl;
    
    
    /*srand(time(0));
    cout << endl;
    cout << "Entrez le seuil de perte c (c doit etre autour de  0  pour une meilleure approximation) : " << endl;
    cin >> c;
    cout << "faites entrer la date future tau du calcul de proba de perte en annees " << endl;
    cin >> tau;
    cout << "faites entrer la valeur de mu (le drift de l'actif sous jacent avant la date tau , si 8% entrez .08) " << endl;
    cin >> mu;
    cout << endl << endl;
    double temps_initial = clock();
    cout << "La proba de perte pour le seuil c de l'estimateur uniforme est egale a : " << EstimUnif(100000, 1000, Prime, S, K, c, Signe, tau, T, mu, r, vol) << endl;
    double temps_final = clock();
    double temps_exe = (temps_final - temps_initial) / CLOCKS_PER_SEC;
    cout << "le temps d'execution est de " << temps_exe << " secondes" << endl;
    cout << endl << endl;
    double temps_initial2 = clock();
    cout << "La proba de perte pour le seuil c de l'estimateur sequentiel est egale a  : " << EstimSequentiel(10000, 2, 10000, Prime, S, K, c, Signe, tau, T, mu, r, vol) << endl;
    double temps_final2 = clock();
    double temps_exe2 = (temps_final2 - temps_initial2) / CLOCKS_PER_SEC;
    cout << "le temps d'execution est de  " << temps_exe2 << " secondes" << endl;
    cout << endl << endl;;
    double temps_initial3 = clock();
    cout << "La proba de perte pour le seuil c de l'estimateur Threshold est égale a : " << EstimThres(1000, -0.1, Prime, S, K, c, Signe, tau, T, mu, r, vol) << endl;
    double temps_final3 = clock();
    double temps_exe3 = (temps_final3 - temps_initial3) / CLOCKS_PER_SEC;
    cout << "Le temps d'execution est de " << temps_exe3 << " secondes" << endl;*/

    int repeter;
    cout << "Voulez vous simuler encore ? Si oui entrez 0" << endl;
    cin >> repeter;
    if (repeter == 0)
    {
        main();
    }
    else
    {
        cout << "Bonne continuation!" << endl;
        return(0);
    }
}

