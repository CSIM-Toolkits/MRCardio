#include "entropy.h"
#include "itkImage.h"
#include <fstream>
//#include "itkImageToHistogramGenerator.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>
#include <string>
#include <math.h>

using namespace std;
entropy::entropy()
{
}

double entropy::Execute(ImageType::Pointer Image, double q_entropy){
    typedef unsigned char       PixelTypen;
    const unsigned int          Dimension = 2;

    typedef itk::Image<PixelTypen, Dimension > ImageTypen;

    typedef itk::Statistics::ScalarImageToHistogramGenerator<
            ImageTypen >   HistogramGeneratorType;

    HistogramGeneratorType::Pointer histogramGenerator =
            HistogramGeneratorType::New();

    histogramGenerator->SetInput(Image);

    histogramGenerator->SetNumberOfBins( 256 );
    histogramGenerator->SetMarginalScale( 10.0 );

    histogramGenerator->SetHistogramMin(  -0.5 );
    histogramGenerator->SetHistogramMax( 255.5 );

    histogramGenerator->Compute();

    typedef HistogramGeneratorType::HistogramType  HistogramType;

    const HistogramType * histogram = histogramGenerator->GetOutput();

    const unsigned int histogramSize = histogram->Size();

    std::cout << "Histogram size " << histogramSize << std::endl;

    t_final = -1;
    int first_bin;
    int last_bin;
    double norm_histo[256];
    double P1[256];
    double P2[256];
    double data[256];
    int l = 0;
    int total = 0;
    int value_out = 0;
    value_out = histogram->GetFrequency(0,0);
    unsigned int bin;
    for( bin=0; bin < histogramSize; bin++ )
    {
        std::cout << "bin = " << bin << " frequency = ";
        std::cout << histogram->GetFrequency( bin, 0 ) << std::endl;
        total += histogram->GetFrequency(bin,0);
        data[l] = histogram->GetFrequency(bin,0);
        l++;
    }
    data[0] = 0;
    total = total - value_out;
    HistogramType::ConstIterator itr = histogram->Begin();
    HistogramType::ConstIterator end = histogram->End();

    for(int i=0;i<256;i++){
        norm_histo[i] = data[i]/total;
    }
    P1[0] = norm_histo[0];
    P2[0] = 1.0 - P1[0];
    for(int i=1;i<256;i++){
        P1[i]= P1[i-1] + norm_histo[i];
        P2[i]= 1.0 - P1[i];

    }

    /* Determine the first non-zero bin */
    first_bin=0;
    for (int i = 0;i<256;i++) {
        if ( !(fabs(P1[i])<2.220446049250313E-16)) {
            first_bin = i;
            break;
        }
    }

    /* Determine the last non-zero bin */
    last_bin=256 - 1;
    for (int i = 256 - 1; i >= first_bin; i-- ) {
        if ( !(fabs(P2[i])<2.220446049250313E-16)) {
            last_bin = i;
            break;
        }
    }
    //q_entropy = (double) SistemaCardio->ui->horizontalScrollBar_2->value()/100.0;
    double max_entrop = -9999999;
    double som_entrop = 0;
    double entrop_a = 0;
    double entrop_b = 0;
    int t;
    if(q_entropy==1.00){
        for(t=first_bin;t<=last_bin;t++){
            entrop_a = 0;
            for(int i=0;i<=t;i++){
                if(data[i]!=0)
                    entrop_a -= (norm_histo[i]/P1[t]) * log (norm_histo[i]/P1[t]);

            }
            entrop_b = 0;
            for(int j=t+1;j<256;j++){
                if(data[j]!=0)
                    entrop_b -= (norm_histo[j]/P2[t]) * log (norm_histo[j]/P2[t]);
            }

            som_entrop = entrop_a + entrop_b;
            if(max_entrop < som_entrop){
                max_entrop = som_entrop;
                t_final = t;
            }

        }
        cout<<first_bin<<endl;
        cout<<last_bin<<endl;
        cout<<"O valor de t: "<<t_final<<endl;
    }

    else{
        for(t=first_bin;t<=last_bin;t++){
            entrop_a = 0;
            for(int i=0;i<=t;i++){
                if(data[i]!=0)
                    //entrop_a += pow((norm_histo[i]/P1[t]),q_entropy);
                    //entrop_a += pow((norm_histo[i]),q_entropy);
                    entrop_a += pow((P1[t]),q_entropy);
            }
            entrop_b = 0;
            for(int j=t+1;j<256;j++){
                if(data[j]!=0)
                    //entrop_b += pow((norm_histo[j]/P2[t]),q_entropy);
                    //entrop_b += pow((norm_histo[j]),q_entropy);
                    entrop_b += pow((P2[t]),q_entropy);
            }
            entrop_a = (1 - entrop_a)/(1 - q_entropy);
            entrop_b = (1 - entrop_b)/(1 - q_entropy);
            som_entrop = entrop_a + entrop_b + ((1 - q_entropy)*(entrop_a*entrop_b));
            if(max_entrop < som_entrop){
                max_entrop = som_entrop;
                t_final = t;
            }

        }
        cout<<first_bin<<endl;
        cout<<last_bin<<endl;
        cout<<"O valor de t - q_entropy: "<<t_final<<endl;
    }
    cout<<"q_entropy: "<<q_entropy<<endl;

    return t_final;

}
