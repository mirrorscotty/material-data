double AA, AB, EaA, EaB; /* Reaction rate parameters */
double A, B, C; /* Viscosity parameters */
double To, Text_hot, Text_cold, t_heat, Tinf; /* Temperature parameters */
double v, L; /* Heat transfer parameters */

/**
 * @brief Reaction rate
 * Calculate the rate of reaction factoring in increase in concentration as a
 * result of the ice crystals forming.
 * @param T Temperature (K)
 * @param c Concentration (kg/L)
 * @return Derivative of concentration with respect to time.
 */
double reaction_rate(double T, double c)
{
        return -AA*exp(-EaA/(R*T))*Xv_water(To)/Xv_water(T)*c;
}

/* Function to return the external temperature of the can. Based on the time,
 * this is either the hot or cold temperature, so that the can is able to be
 * heated and cooled in a single simulation run.
 */
double T_ext(double t)
{
    if(t <= t_heat) {
        return Text_hot;
    } else {
        return Text_cold;
    }
}

/* Return the external temperature if it's not changing. */
double T_inf() {
    return Tinf;
}

/* Return the initial temperature. Used for setting up the temperature inside
 * the can at the start of the simulation.
 */
double T_init(double t)
{
    return To;
}

/* Calculate the viscosity using an Arrhenius-type equation. */
//double mu(double T)
//{
//    return A*pow(10, (B/(T-C)));
//    /* Source: http://en.wikipedia.org/wiki/Viscosity/Viscosity_of_water */
//}

/* Determine the convective heat transfer coefficient */
//double h(double T)
//{
//    double Re, Pr, p, Cp_wat, k;
//   T = T-273.15; /* Convert to Celcius */
//
//    /* Calculate the required data for water */
//   k = 5.7109e-1 + 1.762e-3*T - 6.703e-6*pow(T, 2);
//    if(T >= 0) {
//        Cp_wat = 4.1289 + 9.0864e-5*T - 5.4731e-6*pow(T, 2);
//    } else {
//        Cp_wat = 4.1289 + 5.3062e-3*T - 9.9516e-4*pow(T, 2);
//    }
//    p = 997.18 + 3.1439e-3*T - 3.7574e-3*pow(T, 2);
//
//    T = T + 273.15; /* Convert back to Kelvin to find viscosity */
//    /* Calculate the Reynolds and Prandlt numbers */
//    Re = p*v*L/mu(T);
//    Pr = Cp_wat*mu(T)/k;
//    
//    /* Source: An introduction to Heat and Mass Transfer (Middleman) */
//    return (0.35 + 0.56*pow(Re,0.52))*pow(Pr,0.3)*k/L;
//}
//
/* Test function to spit out a table of data with the x values in one column
 * and the results of a function the other. */
int output_data()
{
    double min, max;
    int points, i;
    const char *file = "out.csv";

    FILE *fp;
    double x, cp_, k_, rho_;

    min = 0;
    max = 0;
    points = 0;
    i = 0;
    fp = NULL;
    x = 0;
    cp_ = 0;
    k_ = 0;
    rho_ = 0;

    fp = fopen(file, "w");
    if(!fp) {
        report_error("Failed to open file for writing.");
    }

    min = 300;
    max = 480;
    points = 100;

    for (i=1; i <= points; i++) {
        x = min+(max-min)*i/points;
        cp_ = Cp(x);
        k_ = k(x);
        rho_ = rho(x);
        fprintf(fp, "%f,%f,%f,%f\n", x, cp_, k_, rho_);
    }

    if(fp) {
        if(fclose(fp) != 0) {
           report_error("Failed to close file.");
           exit(1);           
        }
    }

    return 0;
}

