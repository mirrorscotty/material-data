#include "../unifac/unifac.h"

int main(int argc, char *argv[])
{
    unifac_data *d;
    d = UnifacCreateData();
    UnifacLoadData(d, "unifac1.csv");
    UnifacPrintTable(d);
    return 0;
}

