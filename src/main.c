# include "tsunami.h"





int main(void)
{   

        

        char *meshName = "../data/PacificMedium.txt";
        char *resultBaseName = "output/tsunamiTiny";
        /*         
        tsunamiCompute(10.0,1000,10,meshName,resultBaseName);
        tsunamiAnimate(10.0,1000,10,meshName,resultBaseName);    
        */
        //tsunamiCompute(10.0,1000,10,meshName,resultBaseName);
        //tsunamiAnimate(10.0,1000,10,meshName,resultBaseName);    
        draw();

       

        exit(0);
        return 0;
}