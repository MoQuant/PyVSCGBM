#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <curl/curl.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <time.h>

using namespace boost::property_tree;

// Data Fetchers

std::string key(){
    return "";
}

std::string historical_data(std::string ticker){
    return "https://financialmodelingprep.com/stable/historical-price-eod/light?symbol=" + ticker + "&apikey=" + key();
}

size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* data) {
    size_t totalSize = size * nmemb;
    data->append((char*)contents, totalSize);
    return totalSize;
}

std::string Request(const std::string& url) {
    CURL* curl;
    CURLcode res;
    std::string response;

    curl = curl_easy_init(); 
    if(curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response); 

        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
    }
    return response;
}

std::map<std::string, std::vector<double>> FinancialData(std::string ticker){
    std::map<std::string, std::vector<double>> result;
    std::string response = Request(historical_data(ticker));
    std::stringstream ss(response);
    ptree df;
    read_json(ss, df);
    for(ptree::const_iterator it = df.begin(); it != df.end(); ++it){
        for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
            if(jt->first == "price"){
                result["S"].push_back(jt->second.get_value<double>());
            }
        }
    }
    std::reverse(result["S"].begin(), result["S"].end());
    for(int i = 1; i < result["S"].size(); ++i){
        result["R"].push_back(result["S"][i]/result["S"][i-1] - 1.0);
    }
    return result;
}

// Math Functions

double Weiner(){
    int num = 15;
    double dw = (rand() % (2*num + 1)) - num;
    return dw / 10.0;
}

// Gamma function using Lanczos approximation
double Gamma(double x) {
    // Coefficients for Lanczos approximation
    static const double g = 7.0;
    static const double coef[] = {
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    };

    if (x < 0.5) {
        // Use reflection formula: Gamma(x) = pi / (sin(pi*x) * Gamma(1-x))
        return M_PI / (sin(M_PI * x) * Gamma(1.0 - x));
    }

    x -= 1.0;
    double a = coef[0];
    double t = x + g + 0.5;

    for (int i = 1; i < 9; ++i) {
        a += coef[i] / (x + i);
    }

    return sqrt(2.0 * M_PI) * pow(t, x + 0.5) * exp(-t) * a;
}

double mean(std::vector<double> x){
    double total = 0;
    for(auto & i : x){
        total += i;
    }
    return total / (double) x.size();
}

double stdev(std::vector<double> x){
    double mu = mean(x);
    double total = 0;
    for(auto & i : x){
        total += pow(i - mu, 2);
    }
    total = total / ((double) x.size() - 1);
    return pow(total, 0.5);
}

double Hurst(std::vector<double> ror){
    double mu_ror = mean(ror);
    double sd_ror = stdev(ror);
    double csum = 0;
    std::vector<double> store;

    for(int i = 0; i < ror.size(); ++i){
        csum += ror[i];
        store.push_back(csum - mu_ror);
    }

    std::sort(store.begin(), store.end());

    double Y1 = store[0];
    double Y2 = store[store.size() - 1];

    double hurst = log((Y2 - Y1)/sd_ror)/log(ror.size());

    return hurst;
}

double KT(double t, double H){
    return pow(t, H - 0.5)/Gamma(H + 0.5);
}

double lastElement(std::vector<double> x){
    return x[x.size() - 1];
}

std::map<std::string, double> Quant(std::vector<double> ror, double ds){
    std::map<std::string, double> result;
    int window = 100;
    int lag = 15;
    std::vector<double> hold, store, v0, v1;
    for(int i = window; i < ror.size(); ++i){
        hold = {ror.begin() + i - window, ror.begin()+i};
        store.push_back(pow(stdev(hold), 2));
    }
    double theta = mean(store);
    double sigma = stdev(store);
    v0 = {store.begin() + lag, store.end()};
    v1 = {store.begin(), store.end() - lag};
    double covariance = 0, variance = 0;
    double a1 = mean(v0), b1 = mean(v1);

    for(int i = 0; i < v0.size(); ++i){
        covariance += (v0[i] - a1)*(v1[i] - b1);
        variance += pow(v1[i] - b1, 2);
    }

    double kappa = -log(covariance/variance)/ds;

    result["Kappa"] = kappa;
    result["Theta"] = theta;
    result["Sigma"] = sigma;
    
    return result;
}

int main()
{
    srand(time(NULL));
    
    std::string ticker = "SBUX";
    std::map<std::string, std::vector<double>> data = FinancialData(ticker);

    double H = Hurst(data["R"]); 

    double tf = 1.0 / 365.0;

    int T = 100;
    int IG = 301;
    double dt = tf / T;
    double ds = tf / IG;
    double coef = 0, total = 0;

    std::map<std::string, double> F = Quant(data["R"], ds);

    double SX = data["S"][data["S"].size() - 1];
    double VX = stdev(data["R"]);

    double rho = -0.3;

    double S0 = SX;
    double V0 = VX;
    for(int t = 0; t < T; ++t){
        total = 0;
        for(int s = 0; s < IG; ++s){
            if(s == 0 || s == IG - 1){
                coef = 1.0;
            } else if(s % 2 == 0){
                coef = 2.0;
            } else {
                coef = 4.0;
            }
            V0 = fmax(0.0, V0);
            total += coef * KT(IG - s, H) * F["Kappa"]*(F["Theta"] - V0)*ds + F["Sigma"]*pow(V0, 0.5)*rho*Weiner();
        }
        V0 += total * (1.0/3.0);
        V0 = fmax(V0, 0.0);
        S0 += S0*pow(V0, 0.5)*rho*Weiner();
    }

    std::cout << "Stock Price: " << lastElement(data["S"]) << std::endl;
    std::cout << "Simulated Price: " << S0 << std::endl;  

    return 0;
}


