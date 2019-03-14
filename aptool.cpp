#include "aptool.h"
#include "ui_aptool.h"
#include <iostream>
#include <fstream>
#include <QFileDialog>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QInputDialog>
#include <QDebug>
#include <QPixmap>
#include <QPainter>
#include <QFile>
#include <QTime>
#include <QProgressDialog>
#include <QMessageBox>
#include <qtextstream.h>
#include "s_hull_pro.h"
#include <math.h>
//#include "stdafx.h"
#include <stdlib.h>
#include "interpolation.h"
#include "lwidget.h"

using namespace alglib;

using namespace std;
using namespace cv;

float ReverseFloat( const float inFloat )
{
    float retVal;
    char *floatToConvert = ( char* ) & inFloat;
    char *returnFloat = ( char* ) & retVal;

    // swap the bytes into a temporary buffer
    returnFloat[0] = floatToConvert[3];
    returnFloat[1] = floatToConvert[2];
    returnFloat[2] = floatToConvert[1];
    returnFloat[3] = floatToConvert[0];

    return retVal;
}

static void getHSH(float theta, float phi, float* hweights, int order)
{
    float cosPhi = cos(phi);
    float cosTheta = cos(theta);
    float cosTheta2 = cosTheta * cosTheta;
    hweights[0] = 1/sqrt(2*M_PI);
    hweights[1] = sqrt(6/M_PI)      *  (cosPhi*sqrt(cosTheta-cosTheta2));
    hweights[2] = sqrt(3/(2*M_PI))  *  (-1. + 2.*cosTheta);
    hweights[3] = sqrt(6/M_PI)      *  (sqrt(cosTheta - cosTheta2)*sin(phi));
    if (order >=2)
    {
        hweights[4] = sqrt(30/M_PI)     *  (cos(2.*phi)*(-cosTheta + cosTheta2));
        hweights[5] = sqrt(30/M_PI)     *  (cosPhi*(-1. + 2.*cosTheta)*sqrt(cosTheta - cosTheta2));
        hweights[6] = sqrt(5/(2*M_PI))  *  (1 - 6.*cosTheta + 6.*cosTheta2);
        hweights[7] = sqrt(30/M_PI)     *  ((-1 + 2.*cosTheta)*sqrt(cosTheta - cosTheta2)*sin(phi));
        hweights[8] = sqrt(30/M_PI)     *  ((-cosTheta + cosTheta2)*sin(2.*phi));
    }
    if (order >=3)
    {
        hweights[9]  = 2*sqrt(35/M_PI)	*	(cos(3.0*phi)*pow((cosTheta - cosTheta2), (float)1.5));
        hweights[10] = sqrt(210/M_PI)	*	(cos(2.0*phi)*(-1 + 2*cosTheta)*(-cosTheta + cosTheta2));
        hweights[11] = 2*sqrt(21/M_PI)  *	(cos(phi)*sqrt(cosTheta - cosTheta2)*(1 - 5*cosTheta + 5*cosTheta2));
        hweights[12] = sqrt(7/(2*M_PI)) *	(-1 + 12*cosTheta - 30*cosTheta2 + 20*cosTheta2*cosTheta);
        hweights[13] = 2*sqrt(21/M_PI)  *	(sqrt(cosTheta - cosTheta2)*(1 - 5*cosTheta + 5*cosTheta2)*sin(phi));
        hweights[14] = sqrt(210/M_PI)  *	(-1 + 2*cosTheta)*(-cosTheta + cosTheta2)*sin(2*phi);
        hweights[15] = 2*sqrt(35/M_PI)  *	pow((cosTheta - cosTheta2), (float)1.5)*sin(3*phi);
    }
}

template<class T>struct compare_index
{
    const T base_arr;
    compare_index (const T arr) : base_arr (arr) {}

    bool operator () (int a, int b) const
    {
        return (base_arr[a] < base_arr[b]);
    }
};


void saveRTI_LRGB(QString filename, int W, int H, int type, QString chroma_img){

    char* files[16] = {"h0.bin","h1.bin","h2.bin","h3.bin","h4.bin","h5.bin","h6.bin","h7.bin","h8.bin","h9.bin","h10.bin","h11.bin","h12.bin","h13.bin","h14.bin","h15.bin"};
    float* pBuff = new float[W*H];
    double max[16], min[16];
    float bias[16];
    float scale[16];
    unsigned char c;
    unsigned char* scaledc = new unsigned char[W*H*4];

    cv::Mat image;
    image = cv::imread(chroma_img.toStdString(), CV_LOAD_IMAGE_COLOR);

    for (int i = 0; i < 4; i++){
        ifstream coef (files[i], ios::in | ios::binary);

        coef.read((char*)pBuff,W*H*sizeof(float));
        min[i] = 99999999999;
        max[i] = -99999999999;

        for (int x= 0; x < W*H; x++)
        {
            if(pBuff[x] > max[i]) max[i]=pBuff[x];
            if(pBuff[x] < min[i]) min[i]=pBuff[x];
        }

        scale[i]=(float) (max[i]-min[i]);
        // scale[i]=(float) 1.0+floor((max[i]-min[i]-1)/256);
        bias[i]=(float)(min[i]);// you can change this value
        // delete(pBuff);
        qDebug() <<"minmax "<<min[i]<< ' ' << max[i] <<'\n';
        qDebug() <<"scale "<<float(scale[i])<<'\n';
        qDebug() <<"bias "<<float(bias[i])<<'\n';
        coef.close();
    }

    for (int i = 0; i < 4; i++){
        qDebug(files[i]);
        ifstream coef (files[i], ios::in | ios::binary);

        coef.read((char*)pBuff,W*H*sizeof(float));
        coef.close();
        for (int x = 0;x<W; x++)
            for (int y = 0; y <H; y++)
            {
                scaledc[x+y*W+W*H*i] = (unsigned char)((pBuff[y*W+x] - min[i])/scale[i]);

            }

    }
    qDebug() << "scaledc " << scaledc[100] << " " << pBuff[100] << " " << pBuff[100]*scale[0]/bias[0];
    ofstream outfile; //, testrgb;
    outfile.open(filename.toLatin1(),ios::binary);
    //  testrgb.open("testrgb.txt",ios::app);

    int rtitype, cold,basisterms,basistype, elemsize;

    rtitype=3;
    cold=1; //color dimension
    basisterms=4;
    basistype=2;
    elemsize=1;
    /*********** header *********************************************/
    outfile <<"#HSH1.2"<<"\n";
    outfile <<rtitype<< "\n";
    outfile<<W<<" "<< H<<" "<<cold<<"\n";
    outfile<<basisterms<<" "<< basistype <<" "<< elemsize <<"\n";
    /************** end header ********************/

    outfile.write(( char *)&scale,4*sizeof(float));
    outfile.write(( char *)&bias,4*sizeof(float));

    qDebug() << "sb " << scale[1] << " " << bias[1];
    //       for(int ii=0;ii<9;ii++){
    //           float vrev= ReverseFloat( scale[ii] );
    //            outfile.write(( char *)&vrev,sizeof(float));
    //        }
    //       for(int ii=0;ii<9;ii++){
    //           float vrev= ReverseFloat( bias[ii] );
    //            outfile.write(( char *)&vrev,sizeof(float));
    //        }

    if(!image.data)
    {
        outfile.close();
        return;
    }

    float maxv, vv;
    unsigned char mah;
    for (int y = 0; y <H; y++)
        for (int x = 0;x<W; x++){
            for(int uu=0;uu<3;uu++)
                for (int i = 0; i < 4; i++)
                {
                    c=scaledc[x+y*W+W*H*i] ;//scaledc[i].at<unsigned char>(x,y);

                    outfile.write(( char *)&c,1);
                }



            //             Vec3b val=image.at<Vec3b>(y,x);

            //             maxv=val[0];
            //             for (int i = 1; i < 3; i++)
            //                 if( val[i] > maxv ) maxv=val[i];

            //             for (int i = 0; i < 3; i++){
            //                 mah = (unsigned char) ((float)val[2-i]*255.0/maxv);
            //                 outfile.write(( char *)&mah,1);
            //         }

        }



    //     float maxv, vv;
    //     unsigned char mah;
    //     for (int y = H-1; y >=0; y--)
    //         for (int x = 0;x<W; x++){
    //             Vec3b val=image.at<Vec3b>(y,x);

    //             maxv=val[0];
    //             for (int i = 1; i < 3; i++)
    //                 if(val[i] > maxv) maxv=val[i];

    //             for (int i = 0; i < 3; i++){
    //                 mah =(unsigned char) ((float)val[2-i]*255.0/maxv);
    //                 outfile.write(( char *)&mah,1);
    //             }
    //         }


    delete(scaledc);

    outfile.close();

}


void savePTM_LRGB(QString filename, int W, int H, QString chroma_img){
    //save the result in the LRGB format

    double max[6], min[6];
    int bias[6];
    char* files[6] = {"ptmC1.bin","ptmC2.bin","ptmC3.bin","ptmC4.bin","ptmC5.bin","ptmC6.bin"};
    float scale[6];
    float* pBuff = new float[W*H];
    for (int i = 0; i <= 5; i++){
        ifstream coef (files[i], ios::in | ios::binary);

        coef.read((char*)pBuff,W*H*sizeof(float));
        min[i] = 99999999999;
        max[i] = -99999999999;

        for (int x= 0; x < W*H; x++)
        {
            if(pBuff[x] > max[i]) max[i]=pBuff[x];
            if(pBuff[x] < min[i]) min[i]=pBuff[x];
        }

        scale[i]=(float) 1.0+floor((max[i]-min[i]-1)/256);
        bias[i]=(int)(-min[i]/scale[i]);// you can change this value
        // delete(pBuff);
        qDebug() <<"minmax "<<min[i]<< ' ' << max[i] <<'\n';
        qDebug() <<"scale "<<int(scale[i])<<'\n';
        qDebug() <<"bias "<<int(bias[i])<<'\n';
        coef.close();
    }

    unsigned char* scaledc = new unsigned char[W*H*6];
    unsigned char test;
    //  float testf;
    for (int i = 0; i <= 5; i++){
        qDebug(files[i]);
        ifstream coef (files[i], ios::in | ios::binary);
        // pBuff = new float[W*H];
        coef.read((char*)pBuff,W*H*sizeof(float));
        coef.close();
        /*
        test=(unsigned char)((pBuff[500*W+500]/scale[i])+(float)bias[i]);

        testf = (test-bias[i])*scale[i];
        qDebug() << i << " !!! " << pBuff[500*W+500] << " " << test << " " << testf;
*/
        for (int x = 0;x<W; x++)
            for (int y = 0; y <H; y++)
            {
                scaledc[x+y*W+W*H*i] = (unsigned char)((pBuff[y*W+x]/scale[i])+(float)bias[i]);
            }


        //delete(pBuff);

    }
    delete(pBuff);

    ofstream outfile;
    outfile.open(filename.toLatin1(),ios::binary);
    outfile <<  "PTM_1.2\n";
    outfile <<  "PTM_FORMAT_LRGB\n";
    outfile << W <<"\n";
    outfile << H <<"\n";
    QString num;
    for (int i = 0; i < 6; i++){
        num=QString::number((float)scale[i]);
        outfile << num.toStdString()<<" ";

    } outfile <<'\n';
    for (int i = 0; i < 6; i++){
        num=QString::number((int)bias[i]);
        outfile << num.toStdString()<<' ';
    }outfile <<'\n';

    int offset;
    unsigned char c;

    cv::Mat image;
    image = cv::imread(chroma_img.toStdString(), CV_LOAD_IMAGE_COLOR);

    if(!image.data)
    {
        outfile.close();
        return;
    }
    for (int y = H-1; y >=0; y--)
        for (int x = 0;x<W; x++)
            for (int i = 0; i < 6; i++)
            {
                c=scaledc[x+y*W+W*H*i] ;//scaledc[i].at<unsigned char>(x,y);
                outfile.write(( char *)&c,1);
            }

    float maxv, vv;
    unsigned char mah;
    for (int y = H-1; y >=0; y--)
        for (int x = 0;x<W; x++){
            Vec3b val=image.at<Vec3b>(y,x);

            maxv=val[0];
            for (int i = 1; i < 3; i++)
                if(val[i] > maxv) maxv=val[i];

            for (int i = 0; i < 3; i++){
                mah =(unsigned char) ((float)val[2-i]*255.0/maxv);
                outfile.write(( char *)&mah,1);
            }
        }


    delete(scaledc);

    outfile.close();


}


void savePTM_RGB(QString filename, int W, int H){
    //save the result in the RGB format

    double max[6*3], min[6*3];
    int bias[6];
    char* files[6] = {"ptmC1.bin","ptmC2.bin","ptmC3.bin","ptmC4.bin","ptmC5.bin","ptmC6.bin"};
    float scale[6*3];
    float* pBuff = new float[3*W*H];

    for (int i = 0; i <= 17; i++){
        min[i] = 99999999999;
        max[i] = -99999999999;
    }



    for (int i = 0; i <= 5; i++){
        ifstream coef (files[i], ios::in | ios::binary);
        //for (int k=0; k<3; k++){
        coef.read((char*)pBuff,3*W*H*sizeof(float));


        for (int x= 0; x < 3*W*H; x++)
        {
            if(pBuff[x] > max[i]) max[i]=pBuff[x];
            if(pBuff[x] < min[i]) min[i]=pBuff[x];
            //   if(pBuff[x] > max[i*k+i]) max[i*k+i]=pBuff[x];
            //   if(pBuff[x] < min[i*k+i]) min[i*k+i]=pBuff[x];
        }

        coef.close();
    }

    //  }

    //  for (int k=0; k<3; k++){

    for (int i = 0; i <= 5; i++){
        scale[i]=(float) 1.0+floor((max[i]-min[i]-1)/256);
        bias[i]=(int)(-min[i]/scale[i]);
        //        scale[i*k+i]=(float) 1.0+floor((max[i*k+i]-min[i*k+i]-1)/256);
        //        bias[i*k+i]=(int)(-min[i*k+i]/scale[i*k+i]);// you can change this value
        // delete(pBuff);
        qDebug() <<"minmax "<<min[i]<< ' ' << max[i] <<'\n';
        qDebug() <<"scale "<<int(scale[i])<<'\n';
        qDebug() <<"bias "<<int(bias[i])<<'\n';

    }
    //}

    unsigned char c;

    unsigned char* scaledc = new unsigned char[W*H*6*3];
    unsigned char test;
    float testf;

    for (int i = 0; i <= 5; i++){
        ifstream coef (files[i], ios::in | ios::binary);
        // pBuff = new float[W*H];
        coef.read((char*)pBuff,3*W*H*sizeof(float));
        coef.close();

        //test=(unsigned char)((pBuff[50*W+500]/scale[i])+(float)bias[i]);

        //testf = (test-bias[i])*scale[i];
        // qDebug() << i << " !!! " << pBuff[50*W+50] << " " << test << " " << testf;

        for(int k=0;k<3;k++)
            for (int x = 0;x<W; x++)
                for (int y = 0; y <H; y++)
                {
                    scaledc[x+y*W+i*W*H+k*6*W*H] = (unsigned char)((pBuff[k*W*H+y*W+x]/scale[i])+(float)bias[i]);
                }


        //delete(pBuff);

    }
    delete(pBuff);



    ofstream outfile;
    outfile.open(filename.toLatin1(),ios::binary);
    outfile <<  "PTM_1.2\n";
    outfile <<  "PTM_FORMAT_RGB\n";
    outfile << W <<"\n";
    outfile << H <<"\n";
    QString num;



    for (int i = 0; i < 6; i++){
        num=QString::number((float)scale[i]);
        outfile << num.toStdString()<<" ";

    } outfile <<'\n';


    for (int i = 0; i < 6; i++){
        num=QString::number((int)bias[i]);
        outfile << num.toStdString()<<' ';
    }outfile <<'\n';
    //}

    //    for (int k=0; k<3; k++){

    for (int k=0; k<3; k++){
        for (int y = H-1; y >=0; y--)
            for (int x = 0;x<W; x++)
                for (int i = 0; i < 6; i++)
                {
                    c=scaledc[k*W*H*6+x+y*W+W*H*i] ;//scaledc[i].at<unsigned char>(x,y);
                    outfile.write(( char *)&c,1);
                }
    }
}


apTool::apTool(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::apTool)
{
    ui->setupUi(this);
    QPixmap image=drawLWidget();
    ui->lWidget->setPixmap(image);

    connect(ui->lWidget, SIGNAL(clicked(int,int)), this, SLOT(on_lWidget_clicked(int,int)));
    iw=new ImageView(this);

}

apTool::~apTool()
{
    delete ui;
}

void apTool::on_pushButton_clicked()
{
    QString fileName;
    fileName = QFileDialog::getOpenFileName(this,
                                            tr("Open APH file"));
    ui->fileNameLine->setText(fileName);
}

void apTool::on_processButton_clicked()
{
    // this is the code to fit models over AP data
    // AP types - variable type: 1 - Luminance char , 2 - luminance short , 3 RGB char, 4 RGB short, 0 error
    // direction info - variable dirtype: 1- constant 2 - interpolated

    QString filename = ui->fileNameLine->text();
    QFile file(filename);
    int last = filename.lastIndexOf(QDir::separator());
    QString folder=filename.left(last+1);


    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    QProgressDialog pdialog("Fitting model","",0,100,this);
    pdialog.setWindowModality(Qt::WindowModal);
    pdialog.setCancelButton(0);
    pdialog.setValue(0);
    pdialog.setWindowTitle("Progress Dialog");
    pdialog.show();



    int type=0;
    int dirtype=0;
    int size[2]={0,0};
    int nimg=0;
    Vec3f* dirs;
    float** dircoeffs;
    QString chroma_img;

    // to output binary files with coefficients
    // Tinsae HERE - if more coeff add related output variable
    ofstream outcoef, outcoef1, outcoef2, outcoef3, outcoef4, outcoef5, outcoef6, outcoef7, outcoef8, outcoef9, outcoef10, outcoef11,outcoef12, outcoef13, outcoef14, outcoef15, outcoef16;

    // read AP header

    QTextStream textStream(&file);
    QString line = textStream.readLine();

    if (line.isNull())
        return;

    QStringList parts = line.split(" ");
    if(parts[0] != QString("APA"))
        return;

    while(!(line.isNull())){
        line = textStream.readLine();

        parts = line.split(" ");
        if(parts[0]== "LUMINANCE_TYPE" && parts[1]== "UNSIGNED_CHAR" ){
            qDebug() << "UNSIGNED CHAR";
            type=1;
        }
        if(parts[0]== "LUMINANCE_TYPE" && parts[1]== "UNSIGNED_SHORT" ){
            qDebug() << "UNSIGNED SHORT";
            type=2;
        }

        if(parts[0]== "COLOR_TYPE" && parts[1]== "UNSIGNED_CHAR" ){
            qDebug() << "COLOR UNSIGNED CHAR";
            type=3;
        }
        if(parts[0]== "COLOR_TYPE" && parts[1]== "UNSIGNED_SHORT" ){
            qDebug() << "COLOR UNSIGNED SHORT";
            type=4;
        }


        if(parts[0]== "CHROMA_IMAGE"  ){
            chroma_img = folder+parts[1];
            qDebug() << chroma_img;
        }

        if(parts[0]== "IMAGE_SIZE"  ){
            size[0]=parts[1].toInt();
            size[1]=parts[2].toInt();
            qDebug() << "size " << size[0] << " " << size[1];
        }

        if(parts[0]== "N_IMAGES"){
            nimg = parts[1].toInt();
            qDebug() << "number of images " <<  nimg;

        }
        if(parts[0]== "DIR_TYPE"){
            if(parts[1] == "constant"){
                dirtype=1;
                qDebug() << "Dir type constant";}
            if(parts[1] == "interpolated"){
                dirtype=2;
                qDebug() << "Dir type interpolated";}
        }


        if(parts[0]== "DIRECTION_COEFFICIENTS"){
            if(dirtype==2){
                dirs = new Vec3f[nimg];
                dircoeffs = new float*[nimg];
                for(int i=0;i<nimg;i++)   {
                    line = textStream.readLine();
                    dircoeffs[i]=new float[9];

                    if (line.isNull()){
                        qDebug() << "error";
                        return;
                    }
                    parts = line.split(" ");
                    if(parts.size() < 6){
                        qDebug() << "error";
                        return;
                    }

                    for(int j=0;j<9;j++)
                        dircoeffs[i][j] = parts[j].toFloat();

                    float c_x=size[0]/2;float c_y=size[1]/2;
                    // computing interpolated direction in the image center
                    dirs[i][0]=dircoeffs[i][0]*c_x+dircoeffs[i][1]*c_y+dircoeffs[i][2];
                    dirs[i][1]=dircoeffs[i][3]*c_x+dircoeffs[i][4]*c_y+dircoeffs[i][5];
                    dirs[i][2]=dircoeffs[i][6]*c_x+dircoeffs[i][7]*c_y+dircoeffs[i][8];

                    double nofa=sqrt(  dirs[i][0]*dirs[i][0]+dirs[i][1]*dirs[i][1]+dirs[i][2]*dirs[i][2]);

                   for(int j=0;j<3;j++)
                            dirs[i][j]=dirs[i][j]/nofa;

                }
            }
        }

        if(parts[0] == "LIGHT_DIRECTIONS"){
            if(dirtype==1){
                dirs = new Vec3f[nimg];
                for(int i=0;i<nimg;i++)   {
                    line = textStream.readLine();
                    if (line.isNull()){
                        qDebug() << "error";
                        return;
                    }
                    parts = line.split(" ");
                    if(parts.size() < 3){
                        qDebug() << "error";
                        return;
                    }

                    double nofa=sqrt(parts[0].toFloat()*parts[0].toFloat()+parts[1].toFloat()*parts[1].toFloat()+parts[2].toFloat()*parts[2].toFloat());
                    dirs[i][0]=parts[0].toFloat()/nofa;
                    dirs[i][1]=parts[1].toFloat()/nofa;
                    dirs[i][2]=parts[2].toFloat()/nofa;

                }
            }



        }
    }

    file.close();

    filename.replace(".aph",".apd");
    QFile filed(filename);

    if (!filed.open(QIODevice::ReadOnly))
        return;



    QDataStream in(&filed);

    ui->msgBox->setText("loaded " + filename);


    // now read and process appearance profiles
    vector<unsigned char> apuc;
    vector<unsigned short> apus;

cv:Mat image = cv::imread(chroma_img.toStdString(), CV_LOAD_IMAGE_COLOR);
    //   cv::cvtColor(image,image, cv::COLOR_BGR2RGB);

    Mat test(size[1],size[0],CV_8UC3);
    Mat relight(size[1],size[0],CV_8UC3);
    Mat him(size[1],size[0],CV_8UC3);
    Mat shim(size[1],size[0],CV_8UC3);
    Mat outlim(size[1],size[0],CV_8UC3);
    Mat normals(size[1],size[0],CV_8UC3);
    Mat normals2(size[1],size[0],CV_8UC3);
    normals2=Scalar::all(0);
    Mat albedo(size[1],size[0],CV_8UC1);
    Mat albedog(size[1],size[0],CV_8UC1);
    Mat albedob(size[1],size[0],CV_8UC1);
    std::vector<cv::Mat> albedos(3);

    Vec3f val;
    Vec3f dire;

    unsigned char* valc;
    if(type==1 || type==3) valc = new unsigned char[nimg];

    unsigned short* vals;
    if(type==2 || type==4) vals = new unsigned short[nimg];


    Mat b_l,sol_l;
    Mat L_uv;

    std::vector<Shx> pts, pts2, hull, pt, ptt, hu;
    Shx poi,po;


    // do some tests on light directions sorting and triangulating

    vector<float> elev;
    vector<float> azim;
    vector<float> xx;
    vector<float> yy;

    float dp[nimg][nimg];

    int nhi,nsh, cla;

    pts.clear();

    for (int k=0; k < nimg; k++){


        elev.push_back(asin(dirs[k][2]));
        azim.push_back(atan2(dirs[k][1], dirs[k][0]));

        poi.id = k;
        poi.r = (elev[k]); // pts.txt
        poi.c = ((M_PI+azim[k])/M_PI);
        pts.push_back(poi);

        po.id = k;
        po.r = asin(cos(elev[k]))*sin(azim[k]); // pts.txt
        po.c = asin(cos(elev[k]))*cos(azim[k]);
        pt.push_back(po);
        xx.push_back(po.r);
        yy.push_back(po.c);

        for (int l=0; l<nimg;l++){ // dp difference sqrt(1-(dk*dl)^2)= sin alpha
            dp[k][l] = sqrt(1-(dirs[k][0]*dirs[l][0]+
                    dirs[k][1]*dirs[l][1]+
                    dirs[k][2]*dirs[l][2])*(dirs[k][0]*dirs[l][0]+
                    dirs[k][1]*dirs[l][1]+
                    dirs[k][2]*dirs[l][2]));
        }
    }


    // inds will contain indices in elevation order

    vector<size_t> inds(nimg);
    for (size_t i = 0; i < nimg; ++i) inds[i] = i;
    sort (inds.begin (), inds.end (), compare_index<vector<float> &>(elev));

    for (size_t i = 0; i != inds.size(); ++i) {
        qDebug() << i << " " << inds[i] << " " << elev[inds[i]] << " " << azim[inds[i]];
    }

    // indxx will contain indices in
    vector<size_t> indxx(nimg);
    for (size_t i = 0; i < nimg; ++i) indxx[i] = i;
    sort (indxx.begin (), indxx.end (), compare_index<vector<float> &>(xx));

    vector<int> neigh[nimg];

    std::vector<Triad> triads;

    std::vector<int> outx;

    int nx = de_duplicateX( pts, outx, pts2);
    pts=pts2;

    write_Shx(pts2, "pts.mat");

    int ts = s_hull_pro( pts2, triads);

    write_Triads(triads, "triangles.txt");
    // triangulation in the first space

  /* other triangulation in different space

     std::vector<Triad> triad2;

    std::vector<int> outx2;

    int nx2 = de_duplicateX( pt, outx2, ptt);

    int t2 = s_hull_pro( ptt, triad2);
    // triangulation in the second space
*/
    /* Test saving light points triangulation (works!)*/
    if(0){
        ofstream om;

        om.open ("mesh.off");
        om << "OFF\n";
        om << pts.size() << " " << triads.size() << " " << "\n";

        for(int ii=0;ii<pts.size();ii++){
            om << " " << elev[inds[ii]] << " " << azim[inds[ii]] << " 0 \n";
        }

        for( size_t i = 0; i < triads.size(); i++ )
        {

            om << "3 " << triads[i].a <<  " " << triads[i].b << " " << triads[i].c << endl;

        }
        om.close();
    }
    unsigned char nema[nimg][nimg];
    for (int i=0;i<nimg;i++)
        for (int j=0;j<nimg;j++)
            nema[i][j]=0;

    // adjacency matrix for vertices
    for( size_t i = 0; i < triads.size(); i++ )
    {
        nema[triads[i].a][triads[i].b]=1;
        nema[triads[i].b][triads[i].a]=1;
        nema[triads[i].b][triads[i].c]=1;
        nema[triads[i].c][triads[i].b]=1;
        nema[triads[i].c][triads[i].a]=1;
        nema[triads[i].a][triads[i].c]=1;
    }




    for (int i=0;i<nimg;i++)
        for (int j=0;j<nimg;j++)
            if(nema[i][j]==1)
                neigh[i].push_back(j);

    qDebug() << "nei 0 " << neigh[0][0] << " " << neigh[0][1];


    pts.clear();
    pts2.clear();
    pt.clear();

    // still not really used,
    vector<size_t> il(nimg); //inliers
    vector<size_t> hl(nimg); //highlight
    vector<size_t> sh(nimg); //shadow
    vector<float> shd(nimg); //shadow
    vector<size_t> idx(nimg);


    //   Tinsae HERE - add case for HSH with names like h1.bin, etc.

    if( ui->fitterMenu->currentIndex()==0){
        outcoef.open ("ptmC1.bin", ios::out | ios::binary);
        outcoef1.open ("ptmC2.bin", ios::out | ios::binary);
        outcoef2.open ("ptmC3.bin", ios::out | ios::binary);
        outcoef3.open ("ptmC4.bin", ios::out | ios::binary);
        outcoef4.open ("ptmC5.bin", ios::out | ios::binary);
        outcoef5.open ("ptmC6.bin", ios::out | ios::binary);
    }
    if(ui->fitterMenu->currentIndex()==2){
        outcoef.open ("DrewPtmC1.bin", ios::out | ios::binary);
        outcoef1.open ("DrewPtmC2.bin", ios::out | ios::binary);
        outcoef2.open ("DrewPtmC3.bin", ios::out | ios::binary);
        outcoef3.open ("DrewPtmC4.bin", ios::out | ios::binary);
        outcoef4.open ("DrewPtmC5.bin", ios::out | ios::binary);
        outcoef5.open ("DrewPtmC6.bin", ios::out | ios::binary);
    }
    if(ui->fitterMenu->currentIndex()==1){
        outcoef.open ("albedo.bin", ios::out | ios::binary);
        outcoef1.open ("nx.bin", ios::out | ios::binary);
        outcoef2.open ("ny.bin", ios::out | ios::binary);
        outcoef3.open ("nz.bin", ios::out | ios::binary);
    }

    if(ui->fitterMenu->currentIndex()==3) {
        outcoef.open ("h0.bin", ios::out | ios::binary);
        outcoef1.open ("h1.bin", ios::out | ios::binary);
        outcoef2.open ("h2.bin", ios::out | ios::binary);
        outcoef3.open ("h3.bin", ios::out | ios::binary);
        outcoef4.open ("h4.bin", ios::out | ios::binary);
        outcoef5.open ("h5.bin", ios::out | ios::binary);
        outcoef6.open ("h6.bin", ios::out | ios::binary);
        outcoef7.open ("h7.bin", ios::out | ios::binary);
        outcoef8.open ("h8.bin", ios::out | ios::binary);
        outcoef9.open ("h9.bin", ios::out | ios::binary);
        outcoef10.open ("h10.bin", ios::out | ios::binary);
        outcoef11.open ("h11.bin", ios::out | ios::binary);
        outcoef12.open ("h12.bin", ios::out | ios::binary);
        outcoef13.open ("h13.bin", ios::out | ios::binary);
        outcoef14.open ("h14.bin", ios::out | ios::binary);
        outcoef15.open ("h15.bin", ios::out | ios::binary);
    }

    if(ui->fitterMenu->currentIndex()==4) {return; } // DMD
    if(ui->fitterMenu->currentIndex()==5) {return; }// 3 order ptm

    qDebug() << "type " << type;
    // loop over APA pixel blocks

    if(type <3) {  // if Luminance types, only 1 channel
        for(int j=0;j<size[1];j++)
            for(int i=0;i<size[0];i++)
            {
                pdialog.setValue(100*j/size[1]);
                pdialog.update();
                //            shmask[i][j] = 0;

                if(dirtype==2) // interpolated dirs: estimate pixel-specific direction
                    for(int k=0;k<nimg;k++)
                    {
                        dirs[k][0]=dircoeffs[k][0]*i+dircoeffs[k][1]*j+dircoeffs[k][2];
                        dirs[k][1]=dircoeffs[k][3]*i+dircoeffs[k][4]*j+dircoeffs[k][5];
                        dirs[k][2]=dircoeffs[k][6]*i+dircoeffs[k][7]*j+dircoeffs[k][8];

                        double nofa=sqrt(  dirs[k][0]*dirs[k][0]+dirs[k][1]*dirs[k][1]+dirs[k][2]*dirs[k][2]);

                       for(int uu=0;uu<3;uu++)
                                dirs[k][uu]=dirs[k][uu]/nofa;
                    }


                // 8 bit
                if(type==1){
                    in.readRawData((char*)&valc[0], nimg);
                    vector<unsigned char> vec (&valc[0], &valc[0]+nimg );
                    int medi, mini, maxi;
                    vector<int> shad;
                    vector<int> shad2, tmpv;

                    vector<int> inli;
                    int out;

                    for (size_t p = 0; p != idx.size(); ++p) idx[p] = p;
                    sort (idx.begin (), idx.end (), compare_index<vector<unsigned char> &>(vec));
                    mini = idx[0];
                    maxi = idx[nimg-1];
                    if((vec.size()/2) % 2 == 0)
                        medi = idx[(vec.size())/2 + 1];
                    else
                        medi = idx[(vec.size())/2];


                    nhi=0;
                    nsh=0;
                    for (int p = 0; p < idx.size(); p++) {
                        hl[p]=il[p]=sh[p]=0;
                        shd[p]=0;
                    }

                    // highlights as saturated pixels?
                    for(int k=nimg-1; vec[idx[k]]>0.8*255; k--){nhi++;
                        hl[idx[k]]=1;
                    }

                    for(int k=0;k<nimg;k++){

                        int aa=0;
                        float avgd=0, avgd2=0, avgv=0, avgv2=0;

                        for(int p=0;p<nimg;p++){
                            avgd2 += dp[p][k];
                            avgv2 += sqrt(((float)vec[p]-(float)vec[k])*((float)vec[p]-(float)vec[k]));
                        }


                        //  if(hl[k]==0)
                        for(int p=0; p<neigh[k].size();p++){
                            avgd += dp[neigh[k][p]][k];
                            avgv += sqrt(((float)vec[neigh[k][p]]-(float)vec[k])*((float)vec[neigh[k][p]]-(float)vec[k]));

                            // detect shadows

                            // ombra se vicino è più luminoso in modo doppio a quanto previsto dalla legge del sin elev
                            // edit aggiungo minore della mediana

                            if( vec[k] < vec[medi] &&
                                    (float)vec[neigh[k][p]]*sin(elev[k])/(sin(elev[neigh[k][p]])*(float)vec[k]) > 1.5){
                                sh[k]=1; // c'è almeno un salto: è ombra

                                shd[k]=shd[k]+1; // numero direzioni da cui è in ombra
                                shad.push_back(k);
                                //il[neigh[k][p]]=1;
                                //inli.push_back(neigh[k][p]);
                                // qDebug() << "?? " << vec[inds[k]] << " " << vec[inds[neigh[inds[k]][p]]];
                                aa++;
                            }
                        }
                        avgd2 /= nimg;
                        avgd /=neigh[k].size();
                        avgv2 /= nimg;
                        avgv /=neigh[k].size();

                        if(aa>0)
                            nsh++;
                    }

                    //qDebug() << "sh 1 -" << nsh;
                    out =1;

                    int k=0, ind;

                    for(int k=0;k<3;k++)
                        val[k] = shd[0];

                    test.at<Vec3b>(j,i) = val;


                } // end 8 bit

                if(type==2){   // 16 bit

                    // read pixel data and do some processing. on this skeleton we can
                    // develop fitters, estimate normals, detect shadows and edges, etc.

                    in.readRawData((char*)&vals[0], nimg*2);
                    vector<unsigned short> vec (&vals[0], &vals[0]+nimg );
                    int medi, mini, maxi;
                    vector<int> shad;
                    vector<int> shad2, tmpv;
                    vector<int> inli;

                    int out;

                    for (size_t p = 0; p != idx.size(); ++p) idx[p] = p;
                    sort (idx.begin (), idx.end (), compare_index<vector<unsigned short> &>(vec));

                    mini = idx[0];
                    maxi = idx[nimg-1];
                    if((vec.size()/2) % 2 == 0)
                        medi = idx[(vec.size())/2 + 1];
                    else
                        medi = idx[(vec.size())/2];


                    nhi=0;
                    nsh=0;
                    for (int p = 0; p < idx.size(); p++) {
                        hl[p]=il[p]=sh[p]=0;
                        shd[p]=0;
                        //   qDebug() << vec[idx[p]];
                    }


                    // highlights as saturated pixels?
                    for(int k=nimg-1; vec[idx[k]]>0.8*65535; k--){nhi++;
                        hl[idx[k]]=1;
                    }

#if 0  // only tests exporting files with light positions, ap values etc.
                    // test writing points and values

                    if(i==100 && j==100){

                        QFile tfile("points.txt");
                        if (!tfile.open(QFile::WriteOnly | QFile::Text)) {
                            qDebug() << "error";
                        }
                        QTextStream tstream( &tfile );
                        for(int k=0; k < nimg; k++)
                            tstream << elev[k] << " " << azim[k] << " " << vec[k] << "\n";
                        tfile.close();

                        QFile tfile2("points2.txt");
                        if (!tfile2.open(QFile::WriteOnly | QFile::Text)) {
                            qDebug() << "error";
                        }
                        QTextStream tstream2( &tfile2 );
                        for(int k=0; k < nimg; k++)
                            tstream2 << xx[k] << " " << yy[k] << " " << vec[k] << "\n";
                        tfile2.close();


                        om.open ("meshc.obj");
                        //     om << "COFF\n";
                        //      om << nimg << " " << triads.size() << " " << "\n";

                        float norml=65535;//vec[idx[0]];
                        for(int ii=0;ii<nimg;ii++){
                            om << "v " << elev[inds[ii]] << " " << azim[inds[ii]] << " 0.0 " << vec[inds[ii]]/norml << " " << vec[inds[ii]]/norml << " " << vec[inds[ii]]/norml << " \n";

                        }


                        for( size_t i = 0; i < triads.size(); i++ )
                        {
                            //  om << "3 " << triads[i].a <<  " " << triads[i].b << " " << triads[i].c << endl;
                            om << "f " << triads[i].a+1 <<  " " << triads[i].b+1 << " " << triads[i].c+1 << endl;

                        }
                        om.close();



                        om.open ("meshc2.obj");

                        for(int ii=0;ii<nimg;ii++){
                            om << "v " << xx[indxx[ii]] << " " << yy[indxx[ii]] << " 0.0 " << vec[indxx[ii]]/norml << " " << vec[indxx[ii]]/norml << " " << vec[indxx[ii]]/norml << " \n";

                        }


                        for( size_t i = 0; i < triad2.size(); i++ )
                        {
                            om << "f " << triad2[i].a+1 <<  " " << triad2[i].b+1 << " " << triad2[i].c+1 << endl;

                        }
                        om.close();


                    }
#endif

                    for(int k=0;k<nimg;k++){

                        int aa=0;
                        float avgd=0, avgd2=0, avgv=0, avgv2=0;


                        for(int p=0;p<nimg;p++){
                            avgd2 += dp[p][k];
                            avgv2 += sqrt(((float)vec[p]-(float)vec[k])*((float)vec[p]-(float)vec[k]));
                        }


                        //  if(hl[k]==0)
                        for(int p=0; p<neigh[k].size();p++){
                            avgd += dp[neigh[k][p]][k];
                            avgv += sqrt(((float)vec[neigh[k][p]]-(float)vec[k])*((float)vec[neigh[k][p]]-(float)vec[k]));

                            // detect shadows

                            // ombra se vicino è più luminoso in modo doppio a quanto previsto dalla legge del sin elev
                            // edit aggiungo minore della mediana

                            if( vec[k] < vec[medi] &&
                                    (float)vec[neigh[k][p]]*sin(elev[k])/(sin(elev[neigh[k][p]])*(float)vec[k]) > 3){
                                sh[k]=1; // c'è almeno un salto: è ombra

                                shd[k]=shd[k]+1; // numero direzioni da cui è in ombra
                                shad.push_back(k);
                                //il[neigh[k][p]]=1;
                                //inli.push_back(neigh[k][p]);
                                // qDebug() << "?? " << vec[inds[k]] << " " << vec[inds[neigh[inds[k]][p]]];
                                aa++;
                            }
                        }
                        avgd2 /= nimg;
                        avgd /=neigh[k].size();
                        avgv2 /= nimg;
                        avgv /=neigh[k].size();

                        if(aa>0)
                            nsh++;
                    }

                    //qDebug() << "sh 1 -" << nsh;
                    out =1;

                    int k=0, ind;

                    if(0)
                        while(!shad.empty()) {
                            ind=shad.back();
                            shad.pop_back();
                            for(int p=0; p<neigh[shad[k]].size();p++)
                                if(hl[neigh[shad[k]][p]]==0 && il[neigh[shad[k]][p]]==0 && sh[neigh[shad[k]][p]]==0){
                                    shad2.push_back(neigh[shad[k]][p]);
                                    sh[neigh[shad[k]][p]]=1;
                                    nsh++;
                                }
                            tmpv=shad;
                            shad=shad2;
                            shad2=tmpv;
                        }

                    //qDebug() << "sh 2 -" << nsh;



                    for(int k=0;k<3;k++)
                        val[k] = shd[0];

                    test.at<Vec3b>(j,i) = val;
                }
                // end if depth 16


                int l=0,ninl=0;

                // fit full set
                // Tinsae HERE - add hsh case allocating 9 or 16 values for Luv

                if(ui->robustMenu->currentIndex()==0 ){ // full set: take all the light directons

                    if( ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {
                        L_uv=Mat::zeros(nimg,6,DataType<float>::type);
                        b_l=Mat::zeros(nimg,1,CV_32FC1);
                    }

                    if(ui->fitterMenu->currentIndex()==1){//PS
                        L_uv=Mat::zeros(nimg,3,DataType<float>::type);
                        b_l=Mat::zeros(nimg,1,CV_32FC1);
                    }

                    if(ui->fitterMenu->currentIndex()==3){//hsh 16
                        L_uv=Mat::zeros(nimg,16,DataType<float>::type);
                        b_l=Mat::zeros(nimg,1,CV_32FC1);
                    }

                    if(ui->fitterMenu->currentIndex()==0) {// standard PTM

                        for (int k=0; k<nimg;k++){
                            L_uv.at<float>(k,0)=pow( dirs[k][0] ,2);
                            L_uv.at<float>(k,1)=pow( dirs[k][1],2);
                            L_uv.at<float>(k,2)=(dirs[k][0] * dirs[k][1]);
                            L_uv.at<float>(k,3)=dirs[k][0];
                            L_uv.at<float>(k,4)=dirs[k][1];
                            L_uv.at<float>(k,5)=1;
                            if(type==2)
                                b_l.at<float>(k)= vals[k]/256.0;
                            else if(type==1)
                                b_l.at<float>(k)= (float)valc[k];

                        }
                    }
                    if(ui->fitterMenu->currentIndex()==2) {// Drews PTM

                        for (int k=0; k<nimg;k++){
                            L_uv.at<float>(k,0)=dirs[k][0];
                            L_uv.at<float>(k,1)=dirs[k][1];
                            L_uv.at<float>(k,2)=dirs[k][2];
                            L_uv.at<float>(k,3)=pow( dirs[k][0] ,2);
                            L_uv.at<float>(k,4)=(dirs[k][0] * dirs[k][1]);
                            L_uv.at<float>(k,5)=1;

                            if(type==2)
                                b_l.at<float>(k)= vals[k]/256.0;
                            else if(type==1)
                                b_l.at<float>(k)= (float)valc[k];
                        }
                    }
                    if(ui->fitterMenu->currentIndex()==1) {// PS

                        for (int k=0; k<nimg;k++){
                            L_uv.at<float>(k,0)=dirs[k][0];
                            L_uv.at<float>(k,1)=dirs[k][1];
                            L_uv.at<float>(k,2)=dirs[k][2];

                            if(type==2)
                                b_l.at<float>(k)= vals[k]/256.0;
                            else if(type==1)
                                b_l.at<float>(k)= (float)valc[k];
                        }
                    }

                    if(ui->fitterMenu->currentIndex()==3) {// HSH
                        for (int k=0; k<nimg;k++){
                            float theta=acos(sqrt(1-pow( dirs[k][0] ,2) - pow( dirs[k][1] ,2) )); //angle theta
                            float phi=atan2(dirs[k][0],dirs[k][1]);
                            float hweights[16];
                            float chweights[16];

                            getHSH(theta, phi, hweights, 3);

                            for(int t=0;t<16;t++){
                                L_uv.at<float>(k,t)=hweights[t];
                                //qDebug() << chweights[t];
                            }
                            if(type==2)
                                b_l.at<float>(k)= vals[k]/256.0;
                            else if(type==1)
                                b_l.at<float>(k)= (float)valc[k];
                        }
                    }


                }
                else if( ui->robustMenu->currentIndex()==1 ){  // Trimmed

                    // trimmed fit
                    ninl=nimg-8;
                    //   for (int k=0; k<nimg;k++)
                    //     if(shd[k]==0) ninl=ninl+1;

                    if( ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {
                        L_uv=Mat::zeros(ninl,6,DataType<float>::type);
                        b_l=Mat::zeros(ninl,1,CV_32FC1);
                    }

                    if(ui->fitterMenu->currentIndex()==1){//PS
                        L_uv=Mat::zeros(ninl,3,DataType<float>::type);
                        b_l=Mat::zeros(ninl,1,CV_32FC1);
                    }

                    if(ui->fitterMenu->currentIndex()==3){//hsh 16
                        L_uv=Mat::zeros(ninl,16,DataType<float>::type);
                        b_l=Mat::zeros(ninl,1,CV_32FC1);
                    }

                    if(ui->fitterMenu->currentIndex()==0) {// standard PTM
                        int k=0;
                        for (int p=0; p<nimg;p++)
                            if(idx[p] >3 && idx[p] < nimg-4 )
                            {
                                L_uv.at<float>(k,0)=pow( dirs[idx[p]][0] ,2);
                                L_uv.at<float>(k,1)=pow( dirs[idx[p]][1],2);
                                L_uv.at<float>(k,2)=(dirs[idx[p]][0] * dirs[idx[p]][1]);
                                L_uv.at<float>(k,3)=dirs[idx[p]][0];
                                L_uv.at<float>(k,4)=dirs[idx[p]][1];
                                L_uv.at<float>(k,5)=1;
                                if(type==2)
                                    b_l.at<float>(k)= vals[idx[p]]/256.0;
                                else if(type==1)
                                    b_l.at<float>(k)= (float)valc[idx[p]];
                                k=k+1;

                            }
                    }
                    if(ui->fitterMenu->currentIndex()==2) {// Drews PTM

                        int k=0;
                        for (int p=0; p<nimg;p++)
                            if(idx[p] >3 && idx[p] < nimg-4 ){
                                L_uv.at<float>(k,0)=dirs[idx[p]][0];
                                L_uv.at<float>(k,1)=dirs[idx[p]][1];
                                L_uv.at<float>(k,2)=dirs[idx[p]][2];
                                L_uv.at<float>(k,3)=pow( dirs[idx[p]][0] ,2);
                                L_uv.at<float>(k,4)=(dirs[idx[p]][0] * dirs[idx[p]][1]);
                                L_uv.at<float>(k,5)=1;

                                if(type==2)
                                    b_l.at<float>(k)= vals[idx[p]]/256.0;
                                else if(type==1)
                                    b_l.at<float>(k)= (float)valc[idx[p]];
                                k=k+1;
                            }
                    }
                    if(ui->fitterMenu->currentIndex()==1) {// PS

                        int k=0;
                        for (int p=0; p<nimg;p++)
                            if(idx[p] >3 && idx[p] < nimg-4 ){
                                L_uv.at<float>(k,0)=dirs[idx[p]][0];
                                L_uv.at<float>(k,1)=dirs[idx[p]][1];
                                L_uv.at<float>(k,2)=dirs[idx[p]][2];

                                if(type==2)
                                    b_l.at<float>(k)= vals[idx[p]]/256.0;
                                else if(type==1)
                                    b_l.at<float>(k)= (float)valc[idx[p]];
                                k=k+1;
                            }
                    }
                    //
                    if(ui->fitterMenu->currentIndex()==3) {// HSH
                        for (int k=0; k<nimg-8;k++){
                            float theta=acos(sqrt(1-pow( dirs[idx[k]][0] ,2) - pow( dirs[idx[k]][1] ,2) )); //angle theta
                            float phi=atan2(dirs[idx[k]][0],dirs[idx[k]][1]);
                            float hweights[16];
                            float chweights[16];

                            getHSH(theta, phi, hweights, 3);

                            for(int t=0;t<16;t++){
                                L_uv.at<float>(k,t)=hweights[t];
                            }
                            if(type==2)
                                b_l.at<float>(k)= vals[idx[k]]/256.0;
                            else if(type==1)
                                b_l.at<float>(k)= (float)valc[idx[k]];
                        }
                    }
                    if(ui->fitterMenu->currentIndex()==4) {}// DMD
                    if(ui->fitterMenu->currentIndex()==5) {}// 3 order ptm
                }

                else if(ui->robustMenu->currentIndex()==3 ){ // TEST: removed vs relighted

                    if( ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {
                        L_uv=Mat::zeros(nimg-1,6,DataType<float>::type);
                        b_l=Mat::zeros(nimg-1,1,CV_32FC1);
                    }

                    if(ui->fitterMenu->currentIndex()==1){//PS
                        L_uv=Mat::zeros(nimg-1,3,DataType<float>::type);
                        b_l=Mat::zeros(nimg-1,1,CV_32FC1);
                    }

                    if(ui->fitterMenu->currentIndex()==3){//hsh
                        L_uv=Mat::zeros(nimg-1,16,DataType<float>::type);
                        b_l=Mat::zeros(nimg-1,1,CV_32FC1);
                    }

                    if(ui->fitterMenu->currentIndex()==0) {// standard PTM
                        int k=0;
                        for (int kk=0; kk<nimg;kk++){
                            if(kk!=ui->spinBox->value()){
                                L_uv.at<float>(k,0)=pow( dirs[kk][0] ,2);
                                L_uv.at<float>(k,1)=pow( dirs[kk][1],2);
                                L_uv.at<float>(k,2)=(dirs[kk][0] * dirs[kk][1]);
                                L_uv.at<float>(k,3)=dirs[kk][0];
                                L_uv.at<float>(k,4)=dirs[kk][1];
                                L_uv.at<float>(k,5)=1;
                                if(type==2)
                                    b_l.at<float>(k)= vals[kk]/256.0;
                                else if(type==1)
                                    b_l.at<float>(k)= (float)valc[kk];
                                k++;
                            }
                        }
                    }
                    if(ui->fitterMenu->currentIndex()==2) {// Drews PTM

                        int k=0;
                        for (int kk=0; kk<nimg;kk++)
                            if(kk!=ui->spinBox->value()){

                                L_uv.at<float>(k,0)=dirs[kk][0];
                                L_uv.at<float>(k,1)=dirs[kk][1];
                                L_uv.at<float>(k,2)=dirs[kk][2];
                                L_uv.at<float>(k,3)=pow( dirs[kk][0] ,2);
                                L_uv.at<float>(k,4)=(dirs[kk][0] * dirs[kk][1]);
                                L_uv.at<float>(k,5)=1;

                                if(type==2)
                                    b_l.at<float>(k)= vals[k]/256.0;
                                else if(type==1)
                                    b_l.at<float>(k)= (float)valc[k];
                                k++;
                            }
                    }
                    if(ui->fitterMenu->currentIndex()==1) {// PS

                        int k=0;
                        for (int kk=0; kk<nimg;kk++)
                            if(kk!=ui->spinBox->value()){
                                L_uv.at<float>(k,0)=dirs[kk][0];
                                L_uv.at<float>(k,1)=dirs[kk][1];
                                L_uv.at<float>(k,2)=dirs[kk][2];

                                if(type==2)
                                    b_l.at<float>(k)= vals[kk]/256.0;
                                else if(type==1)
                                    b_l.at<float>(k)= (float)valc[kk];
                                k++;
                            }
                    }

                    if(ui->fitterMenu->currentIndex()==3) {// HSH
                        //
                        int k=0;
                        for (int kk=0; kk<nimg; kk++)
                            if(kk!=ui->spinBox->value()){
                                float theta=acos(sqrt(1-pow( dirs[kk][0] ,2) - pow( dirs[kk][1] ,2) )); //angle theta
                                float phi=atan2(dirs[kk][0],dirs[kk][1]);
                                float hweights[16];
                                float chweights[16];

                                getHSH(theta, phi, hweights, 3);

                                for(int t=0;t<16;t++){
                                    L_uv.at<float>(k,t)=hweights[t];
                                    //qDebug() << chweights[t];
                                }
                                if(type==2)
                                    b_l.at<float>(k)= vals[kk]/256.0;
                                else if(type==1)
                                    b_l.at<float>(k)= (float)valc[kk];
                                k++;
                            }
                    }


                }

                solve(L_uv, b_l, sol_l, DECOMP_SVD);
                L_uv.release();
                b_l.release();

                if(ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {// PTM
                    for (int k=0;k<6;k++){
                        if (cvIsNaN(sol_l.at<float>(k))==1){
                            sol_l.at<float>(k)=0;
                        }

                    }

                    // save coeffs
                                            outcoef.write((char*)&sol_l.at<float>(0),sizeof(float));
                                            outcoef1.write((char*)&sol_l.at<float>(1),sizeof(float));
                                            outcoef2.write((char*)&sol_l.at<float>(2),sizeof(float));
                                            outcoef3.write((char*)&sol_l.at<float>(3),sizeof(float));
                                            outcoef4.write((char*)&sol_l.at<float>(4),sizeof(float));
                                            outcoef5.write((char*)&sol_l.at<float>(5),sizeof(float));
                }

                if(ui->fitterMenu->currentIndex()==2){ // DREW
                     float nf = sqrt(sol_l.at<float>(0)*sol_l.at<float>(0)+sol_l.at<float>(1)*sol_l.at<float>(1)+sol_l.at<float>(2)*sol_l.at<float>(2));

                     for(int k=0;k<3;k++)
                         val[2-k] = 255*0.5*(1+sol_l.at<float>(k)/nf);

                     normals.at<Vec3b>(j,i) = val;
                     albedo.at<unsigned char>(j,i) = 0.6*nf;
                     //qDebug() << val[2] << " ----- " << albedo.at<unsigned char>(j,i);
                 }


                if(ui->fitterMenu->currentIndex()==0 ) {// PTM
                    //relight
                    double reva;
                    int kk= ui->spinBox->value();
                    if(ui->robustMenu->currentIndex()==3){
                        reva = sol_l.at<float>(0)*pow( dirs[kk][0] ,2)+
                                sol_l.at<float>(1)*pow( dirs[kk][1],2)+
                                sol_l.at<float>(2)*(dirs[kk][0] * dirs[kk][1])+
                                sol_l.at<float>(3)*dirs[kk][0]+
                                sol_l.at<float>(4)*dirs[kk][1]+
                                sol_l.at<float>(5);

                        Vec3b col=image.at<Vec3b>(j,i);
                        float nc = col[0]+col[1]+col[2];
                        for(int k=0;k<3;k++){
                            val[k] = std::max(0.0,std::min(255.0,3*col[k]*reva/(nc)));
                        }
                        test.at<Vec3b>(j,i)= val;
                    }
                    else {


                        // estimate relighted
                        double reva;

                        float lx = ui->lxSpinBox->value();
                        float ly = ui->lySpinBox->value();
                        float lz = sqrt(1-lx*lx-ly*ly);
                            reva = sol_l.at<float>(0)*pow( lx ,2)+
                                    sol_l.at<float>(1)*pow( ly,2)+
                                    sol_l.at<float>(2)*(lx * ly)+
                                    sol_l.at<float>(3)*lx+
                                    sol_l.at<float>(4)*ly+
                                    sol_l.at<float>(5);

                            Vec3b col=image.at<Vec3b>(j,i);
                            float nc = col[0]+col[1]+col[2];
                            for(int k=0;k<3;k++){
                                val[k] = std::max(0.0,std::min(255.0,3*col[2-k]*reva/(nc)));
                            }

                        relight.at<Vec3b>(j,i)= val;


                        albedo.at<unsigned char>(j,i) = (unsigned char)sol_l.at<float>(5);

                        float speck=0, shy=0, she=0, stand=0;
                        float diffe,thr;
                        for (int k=0; k<nimg;k++){
                            if(type==2)
                                diffe= vals[k]/256.0;
                            else
                                diffe = valc[k];

                            she = she+4*sh[k];

                            // 12

                            thr= 12;


                            if(ui->fitterMenu->currentIndex()==0)
                                diffe = diffe-(sol_l.at<float>(0)*pow( dirs[k][0] ,2)
                                        +sol_l.at<float>(1)*pow( dirs[k][1],2)
                                        +sol_l.at<float>(2)*(dirs[k][0] * dirs[k][1])
                                        +sol_l.at<float>(3)*dirs[k][0]
                                        +sol_l.at<float>(4)*dirs[k][1]
                                        +sol_l.at<float>(5));
                            else
                                diffe = diffe-(sol_l.at<float>(0)*dirs[k][0]
                                        +sol_l.at<float>(1)*dirs[k][1]
                                        +sol_l.at<float>(2)*dirs[k][2]
                                        +sol_l.at<float>(3)*pow( dirs[k][0] ,2)
                                        +sol_l.at<float>(4)*(dirs[k][0] * dirs[k][1])
                                        +sol_l.at<float>(5));                           

                            stand=stand+diffe*diffe;


                            if(diffe>thr)
                                speck=speck+4;
                            if(diffe<-thr)
                                shy=shy+4;

                        }
                        for(int k=0;k<3;k++)
                            val[k]=speck;
                        him.at<Vec3b>(j,i)= val;
                        for(int k=0;k<3;k++)
                            val[k]=shy;
                        shim.at<Vec3b>(j,i)= val;


                        for(int k=0;k<3;k++)
                            val[k]=25*sqrt(stand)/nimg;

                        test.at<Vec3b>(j,i)= val;

                        for(int k=0;k<3;k++)
                            val[k]=speck+shy;

                        outlim.at<Vec3b>(j,i)= val;
                    }


                }
                else if(ui->fitterMenu->currentIndex()==1){ //PS


                    float lx = ui->lxSpinBox->value();
                    float ly = ui->lySpinBox->value();
                    float lz = sqrt(1-lx*lx-ly*ly);

                    float nf = sqrt(sol_l.at<float>(0)*sol_l.at<float>(0)+sol_l.at<float>(1)*sol_l.at<float>(1)+sol_l.at<float>(2)*sol_l.at<float>(2));
                    float coe=nf;
                    outcoef.write((char*)&coe,sizeof(float));
                    float va = (lx*sol_l.at<float>(0)+ly*sol_l.at<float>(1)+lz*sol_l.at<float>(2));

                    coe=sol_l.at<float>(0)/nf;
                    outcoef1.write((char*)&coe,sizeof(float));
                    coe=sol_l.at<float>(1)/nf;
                    outcoef2.write((char*)&coe,sizeof(float));
                    coe=sol_l.at<float>(2)/nf;
                    outcoef3.write((char*)&coe,sizeof(float));

                    /*     outcoef1.write((char*)&sol_l.at<float>(1),sizeof(float));
                    outcoef2.write((char*)&sol_l.at<float>(2),sizeof(float));
                    outcoef3.write((char*)&sol_l.at<float>(2),sizeof(float));*/

                    for(int k=0;k<3;k++)
                        val[2-k] = 255*0.5*(1+sol_l.at<float>(k)/nf);

                    normals.at<Vec3b>(j,i) = val;
                    albedo.at<unsigned char>(j,i) = 0.6*nf;


                    // estimate relighted
                    Vec3b col=image.at<Vec3b>(j,i);
                    float nc = col[0]+col[1]+col[2];
                    for(int k=0;k<3;k++){
                        val[k] = std::max(0.0,std::min(255.0,double(col[2-k]*va)));
                    }

                    relight.at<Vec3b>(j,i)=val ;

                    // estimating features
                    float speck=0, shy=0, she=0, stand=0;
                    float diffe,thr;

                    for (int k=0; k<nimg;k++){
                        if(type==2)
                            diffe= vals[k]/256.0;
                        else
                            diffe = valc[k];

                        she = she+4*sh[k];

                        thr=12;
                        /*thr= 0.1*(sol_l.at<float>(0)*dirs[k][0]
                            +sol_l.at<float>(1)*dirs[k][1]
                            +sol_l.at<float>(2)*dirs[k][2]);*/

                        diffe = diffe-(sol_l.at<float>(0)*dirs[k][0]
                                +sol_l.at<float>(1)*dirs[k][1]
                                +sol_l.at<float>(2)*dirs[k][2]);

                        stand=stand+diffe*diffe;


                        if(diffe>thr)
                            speck=speck+4;
                        if(diffe<-thr)
                            shy=shy+4;

                    }
                    for(int k=0;k<3;k++)
                        val[k]=speck;
                    him.at<Vec3b>(j,i)= val;
                    for(int k=0;k<3;k++)
                        val[k]=shy;
                    shim.at<Vec3b>(j,i)= val;

                    for(int k=0;k<3;k++)
                        val[k]=20*sqrt(stand)/nimg;

                    test.at<Vec3b>(j,i)= val;

                    for(int k=0;k<3;k++)
                        val[k]=speck+shy;

                    outlim.at<Vec3b>(j,i)= val;



                }
                else if(ui->fitterMenu->currentIndex()==3){ //HSH 16

                    //relight
                    double reva;
                    int kk= ui->spinBox->value();

                    float theta=acos(sqrt(1-pow( dirs[kk][0] ,2) - pow( dirs[kk][1] ,2) )); //angle theta
                    float phi=atan2(dirs[kk][0],dirs[kk][1]);
                    float hweights[16];
                    getHSH(theta, phi, hweights, 3);

                    if(ui->robustMenu->currentIndex()==3){
                        reva = 0;
                        for(int p=0; p<16;p++)
                            reva = reva + sol_l.at<float>(p)*hweights[p];

                        Vec3b col=image.at<Vec3b>(j,i);
                        float nc = col[0]+col[1]+col[2];
                        for(int k=0;k<3;k++){
                            val[k] = std::max(0.0,std::min(255.0,3*col[k]*reva/(nc)));

                        }
                        test.at<Vec3b>(j,i)= val;
                    }
                    else { // save coeffs and images
                        // estimate relighted
                        reva = 0;
                        float lx = ui->lxSpinBox->value();
                        float ly = ui->lySpinBox->value();
                        float lz = sqrt(1-lx*lx-ly*ly);
                        float theta=acos(sqrt(lz)); //angle theta
                        float phi=atan2(lx,ly);
                        float hweights[16];
                        getHSH(theta, phi, hweights, 3);

                        for(int p=0; p<16;p++)
                            reva = reva + sol_l.at<float>(p)*hweights[p];

                        Vec3b col=image.at<Vec3b>(j,i);
                        float nc = col[0]+col[1]+col[2];
                        for(int k=0;k<3;k++){
                            val[k] = std::max(0.0,std::min(255.0,3*col[2-k]*reva/(nc)));

                        }
                        relight.at<Vec3b>(j,i)= val;

                        outcoef.write((char*)&sol_l.at<float>(0),sizeof(float));
                        outcoef1.write((char*)&sol_l.at<float>(1),sizeof(float));
                        outcoef2.write((char*)&sol_l.at<float>(2),sizeof(float));
                        outcoef3.write((char*)&sol_l.at<float>(3),sizeof(float));
                        outcoef4.write((char*)&sol_l.at<float>(4),sizeof(float));
                        outcoef5.write((char*)&sol_l.at<float>(5),sizeof(float));
                        outcoef6.write((char*)&sol_l.at<float>(6),sizeof(float));
                        outcoef7.write((char*)&sol_l.at<float>(7),sizeof(float));
                        outcoef8.write((char*)&sol_l.at<float>(8),sizeof(float));
                        outcoef9.write((char*)&sol_l.at<float>(9),sizeof(float));
                        outcoef10.write((char*)&sol_l.at<float>(10),sizeof(float));
                        outcoef11.write((char*)&sol_l.at<float>(11),sizeof(float));
                        outcoef12.write((char*)&sol_l.at<float>(12),sizeof(float));
                        outcoef13.write((char*)&sol_l.at<float>(13),sizeof(float));
                        outcoef14.write((char*)&sol_l.at<float>(14),sizeof(float));
                        outcoef15.write((char*)&sol_l.at<float>(15),sizeof(float));

                    }


                } ////

                else if(ui->fitterMenu->currentIndex()==4){ // HSH 9?

                }

            }
    }
    else{

        // 8 bit RGB - 16 bit RGB
        if(type==3 || type==4){

            for(int cc=0; cc<3;cc++) { // loop over colors
         //       qDebug() << "test";
                for(int j=0;j<size[1];j++){ // pixel loop
                    pdialog.setValue(100*j*cc/(3*size[1]));
                    pdialog.update();
                    for(int i=0;i<size[0];i++)
                    {

                        if(dirtype==2) // interpolated dirs: estimate pixel-specific direction
                            for(int k=0;k<nimg;k++)
                            {
                                dirs[k][0]=dircoeffs[k][0]*i+dircoeffs[k][1]*j+dircoeffs[k][2];
                                dirs[k][1]=dircoeffs[k][3]*i+dircoeffs[k][4]*j+dircoeffs[k][5];
                                dirs[k][2]=dircoeffs[k][6]*i+dircoeffs[k][7]*j+dircoeffs[k][8];
                            }

                        if(type==3){
                        in.readRawData((char*)&valc[0], nimg);
                     //   vector<unsigned char> vec (&valc[0], &valc[0]+nimg );
                        }
                        if(type==4){
                        in.readRawData((char*)&vals[0], nimg*2);
                      //  vector<unsigned short> vec (&vals[0], &vals[0]+nimg );
                        }

                        int l=0,ninl=0;

                        // fit full set
                        if(ui->robustMenu->currentIndex()==0 ){ // full set
                            if( ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {
                                L_uv=Mat::zeros(nimg,6,DataType<float>::type);
                                b_l=Mat::zeros(nimg,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==1){//PS
                                L_uv=Mat::zeros(nimg,3,DataType<float>::type);
                                b_l=Mat::zeros(nimg,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==3){//hsh 16
                                L_uv=Mat::zeros(nimg,16,DataType<float>::type);
                                b_l=Mat::zeros(nimg,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==0) {// standard PTM

                                for (int k=0; k<nimg;k++){
                                    L_uv.at<float>(k,0)=pow( dirs[k][0] ,2);
                                    L_uv.at<float>(k,1)=pow( dirs[k][1],2);
                                    L_uv.at<float>(k,2)=(dirs[k][0] * dirs[k][1]);
                                    L_uv.at<float>(k,3)=dirs[k][0];
                                    L_uv.at<float>(k,4)=dirs[k][1];
                                    L_uv.at<float>(k,5)=1;
                                    if(type==4)
                                        b_l.at<float>(k)= vals[k]/256.0;
                                    else if(type==3)
                                        b_l.at<float>(k)= (float)valc[k];

                                }
                            }
                            if(ui->fitterMenu->currentIndex()==2) {// Drews PTM

                                for (int k=0; k<nimg;k++){
                                    L_uv.at<float>(k,0)=dirs[k][0];
                                    L_uv.at<float>(k,1)=dirs[k][1];
                                    L_uv.at<float>(k,2)=dirs[k][2];
                                    L_uv.at<float>(k,3)=pow( dirs[k][0] ,2);
                                    L_uv.at<float>(k,4)=(dirs[k][0] * dirs[k][1]);
                                    L_uv.at<float>(k,5)=1;

                                    if(type==4)
                                        b_l.at<float>(k)= vals[k]/256.0;
                                    else if(type==3)
                                        b_l.at<float>(k)= (float)valc[k];
                                }
                            }
                            if(ui->fitterMenu->currentIndex()==1) {// PS

                                for (int k=0; k<nimg;k++){
                                    L_uv.at<float>(k,0)=dirs[k][0];
                                    L_uv.at<float>(k,1)=dirs[k][1];
                                    L_uv.at<float>(k,2)=dirs[k][2];

                                    if(type==4)
                                        b_l.at<float>(k)= vals[k]/256.0;
                                    else if(type==3)
                                        b_l.at<float>(k)= (float)valc[k];
                                }
                            }
                            if(ui->fitterMenu->currentIndex()==3) {// HSH
                                for (int k=0; k<nimg;k++){
                                    float theta=acos(sqrt(1-pow( dirs[k][0] ,2) - pow( dirs[k][1] ,2) )); //angle theta
                                    float phi=atan2(dirs[k][0],dirs[k][1]);
                                    float hweights[16];
                                    float chweights[16];

                                    getHSH(theta, phi, hweights, 3);

                                    getHSH(0, 0, chweights, 3);

                                    for(int t=0;t<16;t++){
                                        L_uv.at<float>(k,t)=hweights[t];
                                        //qDebug() << chweights[t];
                                    }
                                    if(type==4)
                                        b_l.at<float>(k)= vals[k]/256.0;
                                    else if(type==3)
                                        b_l.at<float>(k)= (float)valc[k];
                                }
                            }

                        }
                        else if( ui->robustMenu->currentIndex()==1 ){  // trimmed fit

                            if(type==3){
                                vector<unsigned char> vec (&valc[0], &valc[0]+nimg );
                            for (size_t p = 0; p != idx.size(); ++p) idx[p] = p;
                            sort (idx.begin (), idx.end (), compare_index<vector<unsigned char> &>(vec));
                            }
                            else if(type==4){
                                vector<unsigned short> vec (&vals[0], &vals[0]+nimg );
                                for (size_t p = 0; p != idx.size(); ++p) idx[p] = p;
                                sort (idx.begin (), idx.end (), compare_index<vector<unsigned short> &>(vec));
                            }

                            ninl=nimg-8;

                            if( ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {
                                L_uv=Mat::zeros(ninl,6,DataType<float>::type);
                                b_l=Mat::zeros(ninl,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==1){//PS
                                L_uv=Mat::zeros(ninl,3,DataType<float>::type);
                                b_l=Mat::zeros(ninl,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==3){//hsh 16
                                L_uv=Mat::zeros(ninl,16,DataType<float>::type);
                                b_l=Mat::zeros(ninl,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==0) {// standard PTM

                                for (int k=0; k<nimg-8;k++)
                                {
                                    L_uv.at<float>(k,0)=pow( dirs[idx[k+4]][0] ,2);
                                    L_uv.at<float>(k,1)=pow( dirs[idx[k+4]][1],2);
                                    L_uv.at<float>(k,2)=(dirs[idx[k+4]][0] * dirs[idx[k+4]][1]);
                                    L_uv.at<float>(k,3)=dirs[idx[k+4]][0];
                                    L_uv.at<float>(k,4)=dirs[idx[k+4]][1];
                                    L_uv.at<float>(k,5)=1;
                                    if(type==4)
                                        b_l.at<float>(k)= vals[idx[k+4]]/256.0;
                                    else if(type==3)
                                        b_l.at<float>(k)= (float)valc[idx[k+4]];
                                    k=k+1;

                                }
                            }
                            if(ui->fitterMenu->currentIndex()==2) {// Drews PTM

                                for (int k=0; k<nimg-8; k++){
                                    L_uv.at<float>(k,0)=dirs[idx[k+4]][0];
                                    L_uv.at<float>(k,1)=dirs[idx[k+4]][1];
                                    L_uv.at<float>(k,2)=dirs[idx[k+4]][2];
                                    L_uv.at<float>(k,3)=pow( dirs[idx[k+4]][0] ,2);
                                    L_uv.at<float>(k,4)=(dirs[idx[k+4]][0] * dirs[idx[k+4]][1]);
                                    L_uv.at<float>(k,5)=1;

                                    if(type==4)
                                        b_l.at<float>(k)= vals[idx[k+4]]/256.0;
                                    else if(type==3)
                                        b_l.at<float>(k)= (float)valc[idx[k+4]];
                                    k=k+1;
                                }
                            }
                            if(ui->fitterMenu->currentIndex()==1) {// PS

                                for (int k=0; k<nimg-8;k++){
                                    L_uv.at<float>(k,0)=dirs[idx[k+4]][0];
                                    L_uv.at<float>(k,1)=dirs[idx[k+4]][1];
                                    L_uv.at<float>(k,2)=dirs[idx[k+4]][2];

                                    if(type==4)
                                        b_l.at<float>(k)= vals[idx[k+4]]/256.0;
                                    else if(type==3)
                                        b_l.at<float>(k)= (float)valc[idx[k+4]];
                                    k=k+1;
                                }
                            }

                            if(ui->fitterMenu->currentIndex()==3) {// HSH
                                for (int k=0; k<nimg-8;k++){
                                    float theta=acos(sqrt(1-pow( dirs[idx[k]][0] ,2) - pow( dirs[idx[k]][1] ,2) )); //angle theta
                                    float phi=atan2(dirs[idx[k]][0],dirs[idx[k]][1]);
                                    float hweights[16];
                                    float chweights[16];

                                    getHSH(theta, phi, hweights, 3);

                       //             getHSH(0, 0, chweights, 3);

                                    for(int t=0;t<16;t++){
                                        L_uv.at<float>(k,t)=hweights[t];
                                        //qDebug() << chweights[t];
                                    }
                                    if(type==4)
                                        b_l.at<float>(k)= vals[idx[k]]/256.0;
                                    else if(type==3)
                                        b_l.at<float>(k)= (float)valc[idx[k]];
                                }
                            }

                        }
                        else if( ui->robustMenu->currentIndex()== 3) { // TEST
                      
                            if( ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {
                                L_uv=Mat::zeros(nimg-1,6,DataType<float>::type);
                                b_l=Mat::zeros(nimg-1,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==1){//PS
                                L_uv=Mat::zeros(nimg-1,3,DataType<float>::type);
                                b_l=Mat::zeros(nimg-1,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==3){//hsh 16
                                L_uv=Mat::zeros(nimg-1,16,DataType<float>::type);
                                b_l=Mat::zeros(nimg-1,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==0) {// standard PTM

                                int p=0;
                                for (int k=0; k<nimg;k++)
                                    if(k != ui->spinBox->value()){
                                        L_uv.at<float>(p,0)=pow( dirs[k][0] ,2);
                                        L_uv.at<float>(p,1)=pow( dirs[k][1],2);
                                        L_uv.at<float>(p,2)=(dirs[k][0] * dirs[k][1]);
                                        L_uv.at<float>(p,3)=dirs[k][0];
                                        L_uv.at<float>(p,4)=dirs[k][1];
                                        L_uv.at<float>(p,5)=1;
                                        if(type==4)
                                            b_l.at<float>(p)= vals[k]/256.0;
                                        else if(type==3)
                                            b_l.at<float>(p)= (float)valc[k];

                                        p++;

                                    }
                             
                            }
                            if(ui->fitterMenu->currentIndex()==2) {// Drews PTM
                                int k=0;
                                for (int p=0; p<nimg;p++)
                                    if(p != ui->spinBox->value()){
                                        L_uv.at<float>(k,0)=dirs[p][0];
                                        L_uv.at<float>(k,1)=dirs[p][1];
                                        L_uv.at<float>(k,2)=dirs[p][2];
                                        L_uv.at<float>(k,3)=pow( dirs[p][0] ,2);
                                        L_uv.at<float>(k,4)=(dirs[p][0] * dirs[p][1]);
                                        L_uv.at<float>(k,5)=1;

                                        b_l.at<float>(k)= (float)valc[p];
                                        k++;
                                    }
                            }
                            if(ui->fitterMenu->currentIndex()==1) {// PS

                                int k=0;
                                for (int p=0; p<nimg;p++)
                                    if(p != ui->spinBox->value()){
                                        L_uv.at<float>(k,0)=dirs[p][0];
                                        L_uv.at<float>(k,1)=dirs[p][1];
                                        L_uv.at<float>(k,2)=dirs[p][2];

                                        b_l.at<float>(k)= (float)valc[p];
                                        k++;
                                    }
                            }
                            if(ui->fitterMenu->currentIndex()==3) {// HSH
                                //
                                int k=0;
                                for (int kk=0; kk<nimg; kk++)
                                    if(kk!=ui->spinBox->value()){
                                        float theta=acos(sqrt(1-pow( dirs[kk][0],2) - pow( dirs[kk][1],2) )); //angle theta
                                        float phi=atan2(dirs[kk][0],dirs[kk][1]);
                                        float hweights[16];

                                        getHSH(theta, phi, hweights, 3);

                                        for(int t=0;t<16;t++){
                                            L_uv.at<float>(k,t)=hweights[t];
                                        }
                                        if(type==4)
                                            b_l.at<float>(k)= vals[kk]/256.0;
                                        else if(type==3)
                                            b_l.at<float>(k)= (float)valc[kk];
                                        k++;
                                    }
                            }


                        }

                        solve(L_uv, b_l, sol_l, DECOMP_SVD);
                        L_uv.release();
                        b_l.release();

                        if(ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {// PTM or DrewPTM
                            for (int k=0;k<6;k++){
                                if (cvIsNaN(sol_l.at<float>(k))==1){
                                    sol_l.at<float>(k)=0;
                                }

                            }

                            outcoef.write((char*)&sol_l.at<float>(0),sizeof(float));
                            outcoef1.write((char*)&sol_l.at<float>(1),sizeof(float));
                            outcoef2.write((char*)&sol_l.at<float>(2),sizeof(float));
                            outcoef3.write((char*)&sol_l.at<float>(3),sizeof(float));
                            outcoef4.write((char*)&sol_l.at<float>(4),sizeof(float));
                            outcoef5.write((char*)&sol_l.at<float>(5),sizeof(float));

                        }

                        if(ui->fitterMenu->currentIndex()==0) {// PTM

                            double reva;
                            int kk= ui->spinBox->value();
                            if(ui->robustMenu->currentIndex()==3){ //test
                                reva = sol_l.at<float>(0)*pow( dirs[kk][0] ,2)+
                                        sol_l.at<float>(1)*pow( dirs[kk][1],2)+
                                        sol_l.at<float>(2)*(dirs[kk][0] * dirs[kk][1])+
                                        sol_l.at<float>(3)*dirs[kk][0]+
                                        sol_l.at<float>(4)*dirs[kk][1]+
                                        sol_l.at<float>(5);

                                test.at<Vec3b>(j,i)[cc]= std::max(0.0,std::min(255.0,reva));
                                sol_l.release();
                            }
                            else {  // regular PTM estimate and save features

                                float lx = ui->lxSpinBox->value();
                                float ly = ui->lySpinBox->value();
                                float lz = sqrt(1-lx*lx-ly*ly);

                                reva = sol_l.at<float>(0)*pow( lx ,2)+
                                        sol_l.at<float>(1)*pow( ly,2)+
                                        sol_l.at<float>(2)*(lx*ly)+
                                        sol_l.at<float>(3)*lx+
                                        sol_l.at<float>(4)*ly+
                                        sol_l.at<float>(5);

                                relight.at<Vec3b>(j,i)[cc]= std::max(0.0,std::min(255.0,reva));


                                if(cc==0)
                                    albedo.at<unsigned char>(j,i) = (unsigned char)sol_l.at<float>(5);
                                if(cc==1)
                                    albedog.at<unsigned char>(j,i) = (unsigned char)sol_l.at<float>(5);
                                if(cc==2)
                                    albedob.at<unsigned char>(j,i) = (unsigned char)sol_l.at<float>(5);

                                float speck=0, shy=0, she=0, stand=0;
                                float diffe,thr;
                                for (int k=0; k<nimg;k++){
                                    if(type==4)
                                        diffe= vals[k]/256.0;
                                    else
                                        diffe = valc[k];

                                    she = she+4*sh[k];

                                    thr= 12;

                                    if(ui->fitterMenu->currentIndex()==0)
                                        diffe = diffe-(sol_l.at<float>(0)*pow( dirs[k][0] ,2)
                                                +sol_l.at<float>(1)*pow( dirs[k][1],2)
                                                +sol_l.at<float>(2)*(dirs[k][0] * dirs[k][1])
                                                +sol_l.at<float>(3)*dirs[k][0]
                                                +sol_l.at<float>(4)*dirs[k][1]
                                                +sol_l.at<float>(5));
                                    else
                                        diffe = diffe-(sol_l.at<float>(0)*dirs[k][0]
                                                +sol_l.at<float>(1)*dirs[k][1]
                                                +sol_l.at<float>(2)*dirs[k][2]
                                                +sol_l.at<float>(3)*pow( dirs[k][0] ,2)
                                                +sol_l.at<float>(4)*(dirs[k][0] * dirs[k][1])
                                                +sol_l.at<float>(5));

                                    stand=stand+diffe*diffe;

                                    if(diffe>thr)
                                        speck=speck+4;
                                    if(diffe<-thr)
                                        shy=shy+4;
                                }

                                for(int k=0;k<3;k++)
                                    val[k]=speck;

                                him.at<Vec3b>(j,i)= val;

                                for(int k=0;k<3;k++)
                                    val[k]=shy;

                                shim.at<Vec3b>(j,i)= val;


                                for(int k=0;k<3;k++)
                                    val[k]=25*sqrt(stand)/nimg;

                                test.at<Vec3b>(j,i)= val;

                                for(int k=0;k<3;k++)
                                    val[k]=speck+shy;

                                outlim.at<Vec3b>(j,i)= val;
                                 sol_l.release();
                            }
                        }

                        if(ui->fitterMenu->currentIndex()==1){ //PS


                            float lx = ui->lxSpinBox->value();
                            float ly = ui->lySpinBox->value();
                            float lz = sqrt(1-lx*lx-ly*ly);

                            float va = (lx*sol_l.at<float>(0)+ly*sol_l.at<float>(1)+lz*sol_l.at<float>(2));

                            float nf = sqrt(sol_l.at<float>(0)*sol_l.at<float>(0)+sol_l.at<float>(1)*sol_l.at<float>(1)+sol_l.at<float>(2)*sol_l.at<float>(2));
                            float coe;
                            outcoef.write((char*)&nf,sizeof(float));
                            coe=sol_l.at<float>(0)/nf;
                            outcoef.write((char*)&coe,sizeof(float));
                            coe=sol_l.at<float>(1)/nf;
                            outcoef.write((char*)&coe,sizeof(float));
                            coe=sol_l.at<float>(2)/nf;
                            outcoef.write((char*)&coe,sizeof(float));

                            for(int k=0;k<3;k++)
                                val[2-k] = 255*0.5*(1+sol_l.at<float>(k)/nf);

                            normals.at<Vec3b>(j,i) = val;
                            if(cc==0)
                                albedo.at<unsigned char>(j,i) = 0.6*nf;
                            if(cc==1)
                                albedog.at<unsigned char>(j,i) = 0.6*nf;
                            if(cc==2)
                                albedob.at<unsigned char>(j,i) = 0.6*nf;


                            // estimating features
                            float speck=0, shy=0, she=0, stand=0;
                            float diffe,thr;

                            for (int k=0; k<nimg;k++){
                                if(type==4)
                                    diffe= vals[k]/256.0;
                                else
                                    diffe = valc[k];

                                she = she+4*sh[k];

                                thr=12;
                                /*thr= 0.1*(sol_l.at<float>(0)*dirs[k][0]
                                    +sol_l.at<float>(1)*dirs[k][1]
                                    +sol_l.at<float>(2)*dirs[k][2]);*/

                                diffe = diffe-(sol_l.at<float>(0)*dirs[k][0]
                                        +sol_l.at<float>(1)*dirs[k][1]
                                        +sol_l.at<float>(2)*dirs[k][2]);

                                stand=stand+diffe*diffe;


                                if(diffe>thr)
                                    speck=speck+4;
                                if(diffe<-thr)
                                    shy=shy+4;

                            }
                            for(int k=0;k<3;k++)
                                val[k]=speck;
                            him.at<Vec3b>(j,i)= val;
                            for(int k=0;k<3;k++)
                                val[k]=shy;
                            shim.at<Vec3b>(j,i)= val;

                            for(int k=0;k<3;k++)
                                val[k]=20*sqrt(stand)/nimg;

                            test.at<Vec3b>(j,i)= val;

                            for(int k=0;k<3;k++)
                                val[k]=speck+shy;

                            outlim.at<Vec3b>(j,i)= val;



                            relight.at<Vec3b>(j,i)[cc]= std::max(0.0,std::min(255.0,double(va)));

                        }

                        if(ui->fitterMenu->currentIndex()==2){ // DREW PTM features
                            float nf = sqrt(sol_l.at<float>(0)*sol_l.at<float>(0)+sol_l.at<float>(1)*sol_l.at<float>(1)+sol_l.at<float>(2)*sol_l.at<float>(2));

                            for(int k=0;k<3;k++)
                                val[2-k] = 255*0.5*(1+sol_l.at<float>(k)/nf);

                            normals.at<Vec3b>(j,i) = val;
                            if(cc==0)
                                albedo.at<unsigned char>(j,i) = 0.6*nf;
                            if(cc==1)
                                albedog.at<unsigned char>(j,i) = 0.6*nf;
                            if(cc==2)
                                albedob.at<unsigned char>(j,i) = 0.6*nf;
                            //qDebug() << val[2] << " ----- " << albedo.at<unsigned char>(j,i);
                        }

                        if(ui->fitterMenu->currentIndex()==3){ //HSH

                            //relight
                           double reva;
                            int kk= ui->spinBox->value();

                            float theta=acos(sqrt(1-pow( dirs[kk][0] ,2) - pow( dirs[kk][1] ,2) )); //angle theta
                            float phi=atan2(dirs[kk][0],dirs[kk][1]);
                            float hweights[16];
                            getHSH(theta, phi, hweights, 3);

                           // if(i==100 || i==1450)
                             //   qDebug() << hweights[0] << " " << hweights[0] << " " << hweights[1] << " " << hweights[2] << " " << hweights[3] << " " << hweights[4] << " " << hweights[5] << " ";
                                //qDebug() << j << " -- " << theta << " " << phi;

                            if(ui->robustMenu->currentIndex()==3){
                                reva = 0;
                                for(int p=0; p<16;p++)
                                    reva = reva + sol_l.at<float>(p)*hweights[p];

                                test.at<Vec3b>(j,i)[cc]= std::max(0.0,std::min(255.0,reva));


                                 sol_l.release();
                            }
                            else { // save coeffs and images

                                //relight

                                float lx = ui->lxSpinBox->value();
                                float ly = ui->lySpinBox->value();
                                float lz = sqrt(1-lx*lx-ly*ly);

                               double reva;
                                int kk= ui->spinBox->value();

                                float theta=acos(lz); //angle theta
                                float phi=atan2(lx,ly);
                                float hweights[16];
                                getHSH(theta, phi, hweights, 3);

                                    reva = 0;
                                    for(int p=0; p<16;p++)
                                        reva = reva + sol_l.at<float>(p)*hweights[p];

                                    relight.at<Vec3b>(j,i)[cc]= std::max(0.0,std::min(255.0,reva));

                                        for (int k=0;k<16;k++){
                                            if (cvIsNaN(sol_l.at<float>(k))==1){
                                                sol_l.at<float>(k)=0;
                                                qDebug()<< "NAN";
                                            }}



                                outcoef.write((char*)&sol_l.at<float>(0),sizeof(float));
                                outcoef1.write((char*)&sol_l.at<float>(1),sizeof(float));
                                outcoef2.write((char*)&sol_l.at<float>(2),sizeof(float));
                                outcoef3.write((char*)&sol_l.at<float>(3),sizeof(float));
                                outcoef4.write((char*)&sol_l.at<float>(4),sizeof(float));
                                outcoef5.write((char*)&sol_l.at<float>(5),sizeof(float));
                                outcoef6.write((char*)&sol_l.at<float>(6),sizeof(float));
                                outcoef7.write((char*)&sol_l.at<float>(7),sizeof(float));
                                outcoef8.write((char*)&sol_l.at<float>(8),sizeof(float));
                                outcoef9.write((char*)&sol_l.at<float>(9),sizeof(float));
                                outcoef10.write((char*)&sol_l.at<float>(10),sizeof(float));
                                outcoef11.write((char*)&sol_l.at<float>(11),sizeof(float));
                                outcoef12.write((char*)&sol_l.at<float>(12),sizeof(float));
                                outcoef13.write((char*)&sol_l.at<float>(13),sizeof(float));
                                outcoef14.write((char*)&sol_l.at<float>(14),sizeof(float));
                                outcoef15.write((char*)&sol_l.at<float>(15),sizeof(float));

                            }


                        } ////


                    }
                }
                    normals2=normals2+normals/3.0;
            }
            // end loop over colors

            normals=normals2;
            albedos.at(2) = albedo;
            albedos.at(1) = albedog;
            albedos.at(0) = albedob;


        }

        if(type==14){
            // 16 bit RGB

            for(int cc=0; cc<3;cc++) { // loop over colors
                for(int j=0;j<size[1];j++)
                    for(int i=0;i<size[0];i++)
                    {
                        if(dirtype==2) // interpolated dirs: estimate pixel-specific direction
                            for(int k=0;k<nimg;k++)
                            {
                                dirs[k][0]=dircoeffs[k][0]*i+dircoeffs[k][1]*j+dircoeffs[k][2];
                                dirs[k][1]=dircoeffs[k][3]*i+dircoeffs[k][4]*j+dircoeffs[k][5];
                                dirs[k][2]=dircoeffs[k][6]*i+dircoeffs[k][7]*j+dircoeffs[k][8];
                            }

                        in.readRawData((char*)&vals[0], 2*nimg);
                        vector<unsigned short> vec (&vals[0], &vals[0]+nimg );

                        int l=0,ninl=0;

                        // fit full set
                        if(ui->robustMenu->currentIndex()==0 ){ // full set
                            if( ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {
                                L_uv=Mat::zeros(nimg,6,DataType<float>::type);
                                b_l=Mat::zeros(nimg,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==1){//PS
                                L_uv=Mat::zeros(nimg,3,DataType<float>::type);
                                b_l=Mat::zeros(nimg,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==0) {// standard PTM

                                for (int k=0; k<nimg;k++){
                                    L_uv.at<float>(k,0)=pow( dirs[k][0] ,2);
                                    L_uv.at<float>(k,1)=pow( dirs[k][1],2);
                                    L_uv.at<float>(k,2)=(dirs[k][0] * dirs[k][1]);
                                    L_uv.at<float>(k,3)=dirs[k][0];
                                    L_uv.at<float>(k,4)=dirs[k][1];
                                    L_uv.at<float>(k,5)=1;

                                    b_l.at<float>(k)= vals[k]/256.0;


                                }
                            }
                            if(ui->fitterMenu->currentIndex()==2) {// Drews PTM

                                for (int k=0; k<nimg;k++){
                                    L_uv.at<float>(k,0)=dirs[k][0];
                                    L_uv.at<float>(k,1)=dirs[k][1];
                                    L_uv.at<float>(k,2)=dirs[k][2];
                                    L_uv.at<float>(k,3)=pow( dirs[k][0] ,2);
                                    L_uv.at<float>(k,4)=(dirs[k][0] * dirs[k][1]);
                                    L_uv.at<float>(k,5)=1;

                                    b_l.at<float>(k)= vals[k]/256.0;

                                }
                            }
                            if(ui->fitterMenu->currentIndex()==1) {// PS

                                for (int k=0; k<nimg;k++){
                                    L_uv.at<float>(k,0)=dirs[k][0];
                                    L_uv.at<float>(k,1)=dirs[k][1];
                                    L_uv.at<float>(k,2)=dirs[k][2];

                                    if(type==4)
                                        b_l.at<float>(k)= vals[k]/256.0;
                                    else if(type==3)
                                        b_l.at<float>(k)= (float)valc[k];
                                }
                            }

                        } else if( ui->robustMenu->currentIndex()==1 ){  // trimmed fit
                            ninl=nimg-8;
                            //   for (int k=0; k<nimg;k++)
                            //     if(shd[k]==0) ninl=ninl+1;



                            //  idx[0];

                            if( ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {
                                L_uv=Mat::zeros(ninl,6,DataType<float>::type);
                                b_l=Mat::zeros(ninl,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==1){//PS
                                L_uv=Mat::zeros(ninl,3,DataType<float>::type);
                                b_l=Mat::zeros(ninl,1,CV_32FC1);
                            }

                            if(ui->fitterMenu->currentIndex()==0) {// standard PTM
                                int k=0;
                                for (int p=0; p<nimg;p++)
                                    if(idx[p] >3 && idx[p] < nimg-4 )
                                    {
                                        L_uv.at<float>(k,0)=pow( dirs[k][0] ,2);
                                        L_uv.at<float>(k,1)=pow( dirs[k][1],2);
                                        L_uv.at<float>(k,2)=(dirs[k][0] * dirs[k][1]);
                                        L_uv.at<float>(k,3)=dirs[k][0];
                                        L_uv.at<float>(k,4)=dirs[k][1];
                                        L_uv.at<float>(k,5)=1;

                                        b_l.at<float>(k)= vals[k]/256.0;

                                        k=k+1;

                                    }
                            }
                            if(ui->fitterMenu->currentIndex()==2) {// Drews PTM

                                int k=0;
                                for (int p=0; p<nimg;p++)
                                    if(idx[p] >3 && idx[p] < nimg-4 ){
                                        L_uv.at<float>(k,0)=dirs[k][0];
                                        L_uv.at<float>(k,1)=dirs[k][1];
                                        L_uv.at<float>(k,2)=dirs[k][2];
                                        L_uv.at<float>(k,3)=pow( dirs[k][0] ,2);
                                        L_uv.at<float>(k,4)=(dirs[k][0] * dirs[k][1]);
                                        L_uv.at<float>(k,5)=1;

                                        b_l.at<float>(k)= vals[k]/256.0;

                                        k=k+1;
                                    }
                            }
                            if(ui->fitterMenu->currentIndex()==1) {// PS

                                int k=0;
                                for (int p=0; p<nimg;p++)
                                    if(idx[p] >3 && idx[p] < nimg-4 ){
                                        L_uv.at<float>(k,0)=dirs[k][0];
                                        L_uv.at<float>(k,1)=dirs[k][1];
                                        L_uv.at<float>(k,2)=dirs[k][2];


                                        b_l.at<float>(k)= vals[k]/256.0;

                                        k=k+1;
                                    }
                            }


                        }

                        solve(L_uv, b_l, sol_l, DECOMP_SVD);
                        L_uv.release();
                        b_l.release();

                        if(ui->fitterMenu->currentIndex()==0 || ui->fitterMenu->currentIndex()==2) {// PTM
                            for (int k=0;k<6;k++){
                                if (cvIsNaN(sol_l.at<float>(k))==1){
                                    sol_l.at<float>(k)=0;
                                }

                            }

                            outcoef.write((char*)&sol_l.at<float>(0),sizeof(float));
                            outcoef1.write((char*)&sol_l.at<float>(1),sizeof(float));
                            outcoef2.write((char*)&sol_l.at<float>(2),sizeof(float));
                            outcoef3.write((char*)&sol_l.at<float>(3),sizeof(float));
                            outcoef4.write((char*)&sol_l.at<float>(4),sizeof(float));
                            outcoef5.write((char*)&sol_l.at<float>(5),sizeof(float));


                            if(cc==0)
                                albedo.at<unsigned char>(j,i) = (unsigned char)sol_l.at<float>(5);
                            if(cc==1)
                                albedog.at<unsigned char>(j,i) = (unsigned char)sol_l.at<float>(5);
                            if(cc==2)
                                albedob.at<unsigned char>(j,i) = (unsigned char)sol_l.at<float>(5);

                            float speck=0, shy=0, she=0, stand=0;
                            float diffe,thr;
                            for (int k=0; k<nimg;k++){

                                diffe= vals[k]/256.0;

                                she = she+4*sh[k];

                                thr= 12;


                                if(ui->fitterMenu->currentIndex()==0)
                                    diffe = diffe-(sol_l.at<float>(0)*pow( dirs[k][0] ,2)
                                            +sol_l.at<float>(1)*pow( dirs[k][1],2)
                                            +sol_l.at<float>(2)*(dirs[k][0] * dirs[k][1])
                                            +sol_l.at<float>(3)*dirs[k][0]
                                            +sol_l.at<float>(4)*dirs[k][1]
                                            +sol_l.at<float>(5));
                                else
                                    diffe = diffe-(sol_l.at<float>(0)*dirs[k][0]
                                            +sol_l.at<float>(1)*dirs[k][1]
                                            +sol_l.at<float>(2)*dirs[k][2]
                                            +sol_l.at<float>(3)*pow( dirs[k][0] ,2)
                                            +sol_l.at<float>(4)*(dirs[k][0] * dirs[k][1])
                                            +sol_l.at<float>(5));


                                stand=stand+diffe*diffe;


                                if(diffe>thr)
                                    speck=speck+4;
                                if(diffe<-thr)
                                    shy=shy+4;

                            }
                            for(int k=0;k<3;k++)
                                val[k]=speck;
                            him.at<Vec3b>(j,i)= val;
                            for(int k=0;k<3;k++)
                                val[k]=shy;
                            shim.at<Vec3b>(j,i)= val;


                            for(int k=0;k<3;k++)
                                val[k]=25*sqrt(stand)/nimg;

                            test.at<Vec3b>(j,i)= val;

                            for(int k=0;k<3;k++)
                                val[k]=speck+shy;

                            outlim.at<Vec3b>(j,i)= val;

                        }
                        else if(ui->fitterMenu->currentIndex()==1){

                            outcoef.write((char*)&sol_l.at<float>(0),sizeof(float));
                            outcoef1.write((char*)&sol_l.at<float>(1),sizeof(float));
                            outcoef2.write((char*)&sol_l.at<float>(2),sizeof(float));

                            float nf = sqrt(sol_l.at<float>(0)*sol_l.at<float>(0)+sol_l.at<float>(1)*sol_l.at<float>(1)+sol_l.at<float>(2)*sol_l.at<float>(2));

                            for(int k=0;k<3;k++)
                                val[2-k] = 255*0.5*(1+sol_l.at<float>(k)/nf);

                            normals.at<Vec3b>(j,i) = val;
                            if(cc==0)
                                albedo.at<unsigned char>(j,i) = 0.6*nf;
                            if(cc==1)
                                albedog.at<unsigned char>(j,i) = 0.6*nf;
                            if(cc==2)
                                albedob.at<unsigned char>(j,i) = 0.6*nf;


                            // estimating features
                            float speck=0, shy=0, she=0, stand=0;
                            float diffe,thr;

                            for (int k=0; k<nimg;k++){

                                diffe= vals[k]/256.0;

                                she = she+4*sh[k];

                                thr=12;


                                diffe = diffe-(sol_l.at<float>(0)*dirs[k][0]
                                        +sol_l.at<float>(1)*dirs[k][1]
                                        +sol_l.at<float>(2)*dirs[k][2]);

                                stand=stand+diffe*diffe;


                                if(diffe>thr)
                                    speck=speck+4;
                                if(diffe<-thr)
                                    shy=shy+4;

                            }
                            for(int k=0;k<3;k++)
                                val[k]=speck;
                            him.at<Vec3b>(j,i)= val;
                            for(int k=0;k<3;k++)
                                val[k]=shy;
                            shim.at<Vec3b>(j,i)= val;

                            for(int k=0;k<3;k++)
                                val[k]=20*sqrt(stand)/nimg;

                            test.at<Vec3b>(j,i)= val;

                            for(int k=0;k<3;k++)
                                val[k]=speck+shy;

                            outlim.at<Vec3b>(j,i)= val;

                        }

                        if(ui->fitterMenu->currentIndex()==2){
                            float nf = sqrt(sol_l.at<float>(0)*sol_l.at<float>(0)+sol_l.at<float>(1)*sol_l.at<float>(1)+sol_l.at<float>(2)*sol_l.at<float>(2));

                            for(int k=0;k<3;k++)
                                val[2-k] = 255*0.5*(1+sol_l.at<float>(k)/nf);

                            normals.at<Vec3b>(j,i) = val;
                            if(cc==0)
                                albedo.at<unsigned char>(j,i) = 0.6*nf;
                            if(cc==1)
                                albedog.at<unsigned char>(j,i) = 0.6*nf;
                            if(cc==2)
                                albedob.at<unsigned char>(j,i) = 0.6*nf;
                            //qDebug() << val[2] << " ----- " << albedo.at<unsigned char>(j,i);
                        }


                    }

                normals2=normals2+normals/3.0;
            }
            normals=normals2;
            albedos.at(2) = albedo;
            albedos.at(1) = albedog;
            albedos.at(0) = albedob;
        }
    }

    //
    outcoef.close();
    outcoef1.close();
    outcoef2.close();
    outcoef3.close();

    if(ui->fitterMenu->currentIndex()==2 || ui->fitterMenu->currentIndex()==1){
        outcoef4.close();
        outcoef5.close();
    }

    if(ui->fitterMenu->currentIndex()==4){
        outcoef6.close();
        outcoef7.close();
        outcoef8.close();
        outcoef9.close();
        outcoef10.close();
        outcoef11.close();
        outcoef12.close();
        outcoef13.close();
        outcoef14.close();
        outcoef15.close();
    }

    if(ui->fitterMenu->currentIndex()==3)
        filename.replace(".apd",".rti");
    else
        filename.replace(".apd",".ptm");

    int lastc= filename.lastIndexOf(QDir::separator());
    QString lastname=filename.right(filename.size()-lastc-1);


    //Tinsae HERE - if hsh save related .hsh file

    if(type < 3){

        if(ui->fitterMenu->currentIndex()==3) {// convert saved coeffs to viewable rti
            if(ui->robustMenu->currentIndex()==3){ // TEST

                QString name1 = "relighted_HSH_" + QString::number(ui->spinBox->value()) + ".png";

                 imwrite(name1.toStdString(),test);
             
                cv::cvtColor(test,test, cv::COLOR_BGR2RGB);
                iw->setImage(test);
                iw->show();

                ui->msgBox->setText("saved relighted image");
                ui->msgBox->setText(name1);
            }
            else{

                if(ui->fitViewBox->currentIndex()==5){
                    iw->setImage(relight);
                    iw->show();
                }
            saveRTI_LRGB(lastname,size[0],size[1],1,chroma_img);
            ui->msgBox->setText("saved LRGB .rti file" + lastname);
            }
        }


        if(ui->fitterMenu->currentIndex()==0)
            if(ui->robustMenu->currentIndex()==3){ //TEST

                QString name1 = "relighted_PTM_" + QString::number(ui->spinBox->value()) + ".png";

                 imwrite(name1.toStdString(),test);

                cv::cvtColor(test,test, cv::COLOR_BGR2RGB);
                iw->setImage(test);
                iw->show();

                ui->msgBox->setText("saved relighted image");
            }
            else {// convert saved coeffs to viewable ptm
                savePTM_LRGB(lastname,size[0],size[1],chroma_img);

                ui->msgBox->setText("saved LRGB .ptm file" + lastname);


            }

        if(ui->fitterMenu->currentIndex()==1  || ui->fitterMenu->currentIndex()==2)
            if(ui->robustMenu->currentIndex()==3){

                imwrite("relighted.png",test);
                cv::cvtColor(test,test, cv::COLOR_BGR2RGB);
                iw->setImage(test);
                iw->show();

                ui->msgBox->setText("saved relighted image");
            }
            else {

                imwrite("normals.png",normals);
                imwrite("albedo.png",albedo);
                ui->msgBox->setText("saved normals and albedo images");

                cv::cvtColor(normals,normals, cv::COLOR_BGR2RGB);
                cv::cvtColor(albedo, albedo, cv::COLOR_GRAY2BGR);
                if(ui->fitViewBox->currentIndex()==2){
                    iw->setImage(albedo);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==1){
                    iw->setImage(normals);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==4){
                    iw->setImage(test);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==3){
                    iw->setImage(outlim);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==5){
                    iw->setImage(relight);
                    iw->show();
                }

            }

        if(ui->fitterMenu->currentIndex()==0)
            if(ui->robustMenu->currentIndex()==3){
                iw->setImage(test);
                iw->show();
            }
            else {

                if(ui->fitViewBox->currentIndex()==1){
                    ui->msgBox->setText("Normals not estimated in PTM mode");
                }
                if(ui->fitViewBox->currentIndex()==2){
                    ui->msgBox->setText("Albedo not estimated in PTM mode");
                }
                if(ui->fitViewBox->currentIndex()==4){
                    iw->setImage(test);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==3){
                    iw->setImage(outlim);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==5){
                    iw->setImage(relight);
                    iw->show();
                }

            }



        cv::Mat patch, patchi;
        if(ui->fitterMenu->currentIndex()<3 && ui->robustMenu->currentIndex()!= 3){
            imwrite("him.png",him);
            imwrite("shim.png",shim);
            imwrite("outlim.png",outlim);
            imwrite("residual.png",test);

            for(int i=0; i<him.cols; i+=21)
                for(int j=0; j<him.rows; j+=21){
                    cv::Rect roi(i,j,20, 20);
                    patch = test(roi);
                    patchi = albedo(roi);
                    Scalar m,d;
                    meanStdDev(patch, m, d);
                    //pRoi.setTo(cv::Scalar(blue, green, red));
                    //image_roi2.copyTo( image_roi );
                }
        }

    }
    else{ // COLOR - type>2


        if(ui->fitterMenu->currentIndex()==3) {// convert saved coeffs to viewable rti
            if(ui->robustMenu->currentIndex()==3){ // TEST

                iw->setImage(test);
                iw->show();
                cv::cvtColor(test,test, cv::COLOR_BGR2RGB);
                // imwrite("relighted.png",test);
                  QString name1 = "relighted_HSH_" + QString::number(ui->spinBox->value()) + ".png";

                 imwrite(name1.toStdString(),test);


                ui->msgBox->setText("saved relighted image");
            }
            else{

                if(ui->fitViewBox->currentIndex()==5){
                    iw->setImage(relight);
                    iw->show();
                }
          //  saveRTI_LRGB(lastname,size[0],size[1],1,chroma_img);
          //  ui->msgBox->setText("saved LRGB .rti file" + lastname);
            }
        }



        if(ui->fitterMenu->currentIndex()==0){

            // convert saved coeffs to viewable ptm

            savePTM_RGB(lastname,size[0],size[1]);


            // cv::Mat colorImage;


//            if(ui->fitViewBox->currentIndex()==4){
//                iw->setImage(test);
//                iw->show();
//            }
//            if(ui->fitViewBox->currentIndex()==3){
//                iw->setImage(outlim);
//                iw->show();
//            }

            // cv::merge(albedos, colorImage);
            // imwrite("albedo.png",colorImage);

            //            cv::cvtColor(colorImage, colorImage, cv::COLOR_BGR2RGB);
            //            if(ui->fitViewBox->currentIndex()==2){
            //                iw->setImage(colorImage);
            //                iw->show();
            //            }
        }

        if(ui->fitterMenu->currentIndex()==1  || ui->fitterMenu->currentIndex()==2){//PS o drew
            if(ui->robustMenu->currentIndex()==3){



                iw->setImage(test);
                iw->show();
               imwrite("relighted.png",test);
            }
            else {

                imwrite("normals.png",normals);

                cv::Mat colorImage;
                cv::merge(albedos, colorImage);
                imwrite("albedo.png",colorImage);
                //imwrite("albedo.png",albedo);
                ui->msgBox->setText("saved normals and albedo images");


                cv::cvtColor(normals,normals, cv::COLOR_BGR2RGB);
                //cv::cvtColor(albedo, albedo, cv::COLOR_GRAY2BGR);
                if(ui->fitViewBox->currentIndex()==2){
                    iw->setImage(colorImage);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==1){
                    iw->setImage(normals);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==4){
                    iw->setImage(test);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==3){
                    iw->setImage(outlim);
                    iw->show();
                }


                if(ui->fitViewBox->currentIndex()==5){
                    iw->setImage(relight);
                    iw->show();
                }
                imwrite("him.png",him);
                imwrite("shim.png",shim);
                imwrite("outlim.png",outlim);
                imwrite("residual.png",test);
            }
        }

        if(ui->fitterMenu->currentIndex()==0)
        {
            if(ui->robustMenu->currentIndex()==3){
                iw->setImage(test);
                iw->show();
                ui->msgBox->setText("saved relighted image");
                  QString name1 = "relighted_PTM_" + QString::number(ui->spinBox->value()) + ".png";


                cv::cvtColor(test,test, cv::COLOR_BGR2RGB);
                 imwrite(name1.toStdString(),test);
               // imwrite("relighted.png",test);
            }
            else{

                if(ui->fitViewBox->currentIndex()==1){
                    ui->msgBox->setText("Normals not estimated in PTM mode");
                }
                if(ui->fitViewBox->currentIndex()==2){
                    ui->msgBox->setText("Albedo not estimated in PTM mode");
                }
                if(ui->fitViewBox->currentIndex()==4){
                    iw->setImage(test);
                    iw->show();
                }
                if(ui->fitViewBox->currentIndex()==3){
                    iw->setImage(outlim);
                    iw->show();
                }

                if(ui->fitViewBox->currentIndex()==5){
                    iw->setImage(relight);
                    iw->show();
                }
                imwrite("him.png",him);
                imwrite("shim.png",shim);
                imwrite("outlim.png",outlim);
                imwrite("residual.png",test);
            }
        }



    }

    filed.close();
}

void apTool::on_showButton_clicked()   // Funzione per relighting
{
    iw->white1->setGeometry(QRect(0,0,0,0));
    iw->point1->setGeometry(QRect(0,0,0,0));
    iw->CoordinatesSet1=false;
    iw->white2->setGeometry(QRect(0,0,0,0));
    iw->point2->setGeometry(QRect(0,0,0,0));
    iw->CoordinatesSet2=false;

    QString filename = ui->fileNameLine->text();
    QFile file(filename);
    qDebug() << filename;

    int last= filename.lastIndexOf(QDir::separator());
    QString folder=filename.left(last+1);

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;


    QProgressDialog pdialog("Creating image","",0,100,this);
    pdialog.setWindowModality(Qt::WindowModal);
    pdialog.setCancelButton(0);
    pdialog.setValue(0);
    pdialog.setWindowTitle("Progress Dialog");
    pdialog.show();


    int type=0;
    int dirtype=0;
    int size[2]={0,0};
    int nimg=0;
    Vec3f* dirs;
    float** dircoeffs;
    QString chroma_img;


    // Reading header file

    QTextStream textStream(&file);

    QString line = textStream.readLine();

    if (line.isNull())
        return;

    QStringList parts = line.split(" ");
    if(parts[0] != QString("APA"))
        return;

    while(!(line.isNull())){
        line = textStream.readLine();
        // qDebug() << line;

        parts = line.split(" ");
        if(parts[0]== "LUMINANCE_TYPE" && parts[1]== "UNSIGNED_CHAR" ){
            qDebug() << "UNSIGNED CHAR";
            type=1;
        }
        if(parts[0]== "LUMINANCE_TYPE" && parts[1]== "UNSIGNED_SHORT" ){
            qDebug() << "UNSIGNED SHORT";
            type=2;
        }

        if(parts[0]== "COLOR_TYPE" && parts[1]== "UNSIGNED_CHAR" ){
            qDebug() << "COLOR UNSIGNED CHAR";
            type=3;
        }
        if(parts[0]== "COLOR_TYPE" && parts[1]== "UNSIGNED_SHORT" ){
            qDebug() << "COLOR UNSIGNED SHORT";
            type=4;
        }



        if(parts[0]== "CHROMA_IMAGE"  ){
            //  chroma_img = parts[1];
            chroma_img = folder+parts[1];
            qDebug() << chroma_img;

        }

        if(parts[0]== "IMAGE_SIZE"  ){
            size[0]=parts[1].toInt();
            size[1]=parts[2].toInt();
            qDebug() << "size " << size[0] << " " << size[1];
        }

        if(parts[0]== "N_IMAGES"){
            nimg = parts[1].toInt();
            qDebug() << "number of images " <<  nimg;

        }
        if(parts[0]== "DIR_TYPE"){
            if(parts[1] == "constant"){
                dirtype=1;
                qDebug() << "Dir type constant";}
            if(parts[1] == "interpolated"){
                dirtype=2;
                qDebug() << "Dir type interpolated";}
        }


        if(parts[0]== "DIRECTION_COEFFICIENTS"){
            if(dirtype==2){
                dirs = new Vec3f[nimg];
                dircoeffs = new float*[nimg];
                for(int i=0;i<nimg;i++)   {
                    line = textStream.readLine();
                    dircoeffs[i]=new float[9];

                    if (line.isNull()){
                        qDebug() << "error";
                        return;
                    }
                    parts = line.split(" ");
                    if(parts.size() < 6){
                        qDebug() << "error";
                        return;
                    }

                    for(int j=0;j<9;j++)
                        dircoeffs[i][j] = parts[j].toFloat();

                    //qDebug() << dircoeffs[i][0] << " " << dircoeffs[i][1] << " " << dircoeffs[i][5];
                    float c_x=size[0]/2;float c_y=size[1]/2;
                    dirs[i][0]=dircoeffs[i][0]*c_x+dircoeffs[i][1]*c_y+dircoeffs[i][2];
                    dirs[i][1]=dircoeffs[i][3]*c_x+dircoeffs[i][4]*c_y+dircoeffs[i][5];
                    dirs[i][2]=dircoeffs[i][6]*c_x+dircoeffs[i][7]*c_y+dircoeffs[i][8];

                    double nofa=sqrt(  dirs[i][0]*dirs[i][0]+dirs[i][1]*dirs[i][1]+dirs[i][2]*dirs[i][2]);

                   for(int j=0;j<3;j++)
                            dirs[i][j]=dirs[i][j]/nofa;

                    qDebug() << dirs[i][0] << " " << dirs[i][1] << " "  << dirs[i][2];

                    // computing interpolated direction in the center
                    //
                }
            }
        }


        //if light directions are constant
        if(parts[0] == "LIGHT_DIRECTIONS"){
            if(dirtype==1){
                dirs = new Vec3f[nimg];
                for(int i=0;i<nimg;i++)   {
                    line = textStream.readLine();
                    if (line.isNull()){
                        qDebug() << "error";
                        return;
                    }
                    parts = line.split(" ");
                    if(parts.size() < 3){
                        qDebug() << "error";
                        return;
                    }

                    dirs[i][0]=parts[0].toFloat();
                    dirs[i][1]=parts[1].toFloat();
                    dirs[i][2]=parts[2].toFloat();

                    double nofa=sqrt(  dirs[i][0]*dirs[i][0]+dirs[i][1]*dirs[i][1]+dirs[i][2]*dirs[i][2]);

                   for(int j=0;j<3;j++)
                            dirs[i][j]=dirs[i][j]/nofa;

                }
            }

        }
    }

    file.close();


    // Image for output
    Mat shownImg(size[1],size[0],CV_8UC3,cv::Scalar(0,0,0));

    Mat comparedImg(size[1],size[0],CV_8UC3,cv::Scalar(0,0,0));

    ui->msgBox->setText("Estimating image");

    Vec3f val;
    Vec3f valcomp;
    Vec3f dire;
    vector <double> dist;
    unsigned char* valc = new unsigned char[nimg];
    unsigned short* vals = new unsigned short[nimg];
    unsigned char** colc = new unsigned char*[3];
    for(int k=0; k<3; k++){
        colc[k] = new unsigned char[nimg];
    }

    // this was processing to analyze the light directions
    //estimate triangulations of the lights in the central pixel not used now but will use in new smart algorithms...

    std::vector<Shx> pts, pts2, hull, pt, ptt, hu;
    Shx poi,po;

    vector<float> elev;
    vector<float> azim;
    vector<float> xx;
    vector<float> yy;

    float dp[nimg][nimg];

    int nhi,nsh, cla;

    pts.clear();

    for (int k=0; k < nimg; k++){
        elev.push_back(asin(dirs[k][2]));
        azim.push_back(atan2(dirs[k][1], dirs[k][0]));

        poi.id = k;
        poi.r = (elev[k]); // pts.txt
        poi.c = ((M_PI+azim[k])/M_PI);
        pts.push_back(poi);

        po.id = k;
        po.r = dirs[k][2]*sin(azim[k]); // pts.txt
        po.c = dirs[k][2]*cos(azim[k]);
        pt.push_back(po);
        xx.push_back(po.r);
        yy.push_back(po.c);

        for (int l=0; l<nimg;l++){
            dp[k][l] = sqrt(1-(dirs[k][0]*dirs[l][0]+
                    dirs[k][1]*dirs[l][1]+
                    dirs[k][2]*dirs[l][2])*(dirs[k][0]*dirs[l][0]+
                    dirs[k][1]*dirs[l][1]+
                    dirs[k][2]*dirs[l][2]));
        }
    }


    // inds will contain indices in elevation order

    vector<size_t> inds(nimg);
    for (size_t i = 0; i < nimg; ++i) inds[i] = i;
    sort (inds.begin (), inds.end (), compare_index<vector<float> &>(elev));

    for (size_t i = 0; i != inds.size(); ++i) {
        qDebug() << i << " " << inds[i] << " " << elev[inds[i]] << " " << azim[inds[i]];
    }

    // inds will contain indices in lx order
    vector<size_t> indxx(nimg);
    for (size_t i = 0; i < nimg; ++i) indxx[i] = i;
    sort (indxx.begin (), indxx.end (), compare_index<vector<float> &>(xx));

    vector<int> neigh[nimg];

    std::vector<Triad> triads;

    std::vector<int> outx;

    int nx = de_duplicateX( pts, outx, pts2);
    pts=pts2;

    int ts = s_hull_pro( pts2, triads);

    std::vector<Triad> triad2;

    std::vector<int> outx2;

    int nx2 = de_duplicateX( pt, outx2, ptt);

    int t2 = s_hull_pro( ptt, triad2);

    // End of triangulations

    // other stuff to ignore: estimate the neighbors of each light direction (alwasy on the central pixel)
    unsigned char nema[nimg][nimg];
    for (int i=0;i<nimg;i++)
        for (int j=0;j<nimg;j++)
            nema[i][j]=0;

    for( size_t i = 0; i < triads.size(); i++ )
    {
        nema[triads[i].a][triads[i].b]=1;
        nema[triads[i].b][triads[i].a]=1;
        nema[triads[i].b][triads[i].c]=1;
        nema[triads[i].c][triads[i].b]=1;
        nema[triads[i].c][triads[i].a]=1;
        nema[triads[i].a][triads[i].c]=1;

    }

    for (int i=0;i<nimg;i++)
        for (int j=0;j<nimg;j++)
            if(nema[i][j]==1)
                neigh[i].push_back(j);
    pts.clear();
    pts2.clear();
    pt.clear();


    // Da qui la parte che interessa

    // read the chromaticity image in a opencv image


    cv::Mat image;

    if(type < 3){
        image = cv::imread(chroma_img.toStdString(), CV_LOAD_IMAGE_COLOR);
        cv::cvtColor(image,image, cv::COLOR_BGR2RGB);
    }

    iw->setImage(shownImg);
    iw->show();


    // now read the data file
    filename.replace(".aph",".apd");
    QFile filed(filename);

    if (!filed.open(QIODevice::ReadOnly))
        return;

    QDataStream in(&filed);

    bool flagInterp = false;
    rbfmodel modelInterpolation;
    rbfreport repInterpolation;
    real_2d_array arrayInterpolation;
    vector<size_t> iv;

    double lx = ui->lxSpinBox->value();
    double ly = ui->lySpinBox->value();
    double lz = sqrt(1-lx*lx-ly*ly);

    if(ui->viewBox->currentIndex()==2){


        for(int k=0;k<nimg;k++){
            double dv=(lx-dirs[k][0])*(lx-dirs[k][0])+(ly-dirs[k][1])*(ly-dirs[k][1])+(lz-dirs[k][2])*(lz-dirs[k][2]);
            dist.push_back(dv);
        }

        iv.resize(nimg);
        for (size_t p = 0; p != iv.size(); ++p) iv[p] = p;
        sort (iv.begin (), iv.end (), compare_index<vector<double> &>(dist));
        dist.clear();

    }

    if(ui->viewBox->currentIndex()==4){

        int hh=ui->spinBox->value();
        lx=dirs[ui->spinBox->value()][0];
        ly=dirs[ui->spinBox->value()][1];
        lz = sqrt(1-lx*lx-ly*ly);

        dist.clear();
        for(int k = 0; k < nimg; k++){
            double dv=(lx-dirs[k][0])*(lx-dirs[k][0])+(ly-dirs[k][1])*(ly-dirs[k][1])+(lz-dirs[k][2])*(lz-dirs[k][2]);
            dist.push_back(dv);
        }

        iv.resize(nimg);
        for (size_t p = 0; p != iv.size(); ++p) iv[p] = p;
        sort (iv.begin (), iv.end (), compare_index<vector<double> &>(dist));
        dist.clear();

    }


    // Part to look at: loop over APA pixel blocks

    if(type <3)
        for(int j=0;j<size[1];j++)
            for(int i=0;i<size[0];i++)
            {

                if(dirtype==2) // interpolated dirs: estimate pixel-specific direction
                    for(int k=0;k<nimg;k++)
                    {
                        dirs[k][0]=dircoeffs[k][0]*i+dircoeffs[k][1]*j+dircoeffs[k][2];
                        dirs[k][1]=dircoeffs[k][3]*i+dircoeffs[k][4]*j+dircoeffs[k][5];
                        dirs[k][2]=dircoeffs[k][6]*i+dircoeffs[k][7]*j+dircoeffs[k][8];

                        if(i==0 && j==0)
                            qDebug() << " 0,0 - " << dirs[k][0] << " " << dirs[k][1] << " "  << dirs[k][2];
                        if(i==size[0]-1 && j==size[1]-1)
                            qDebug() << " end  - " << dirs[k][0] << " " << dirs[k][1] << " "  << dirs[k][2];
                    }

                // 8 bit LRGB
                if(type==1){
                    in.readRawData((char*)&valc[0], nimg);
                    vector<unsigned char> vec (&valc[0], &valc[0]+nimg );
                    int medi, mini, maxi;
                    vector<int> shad;
                    vector<int> shad2, tmpv;
                    // 8 bit version,

                    // read pixel data and do some processing. on this skeleton we can
                    // develop fitters, estimate normals, detect shadows and edges, etc.

                    float avg;
                    vector<size_t> idx(vec.size());

                    /*  Vec3f cc; cc[0] = 0.2126; cc[1] = 0.7152; cc[2] = 0.0722;*/

                    if(ui->viewBox->currentIndex()==3){ // max
                        for (size_t p = 0; p != idx.size(); ++p) {idx[p] = p; }
                        sort (idx.begin (), idx.end (), compare_index<vector<unsigned char> &>(vec));
                        mini = idx[0];
                        maxi = idx[nimg-1];

                        Vec3b col=image.at<Vec3b>(j,i);
                        float nc = col[0]+col[1]+col[2];
                        for(int k=0;k<3;k++){
                            val[k] = 3*col[k]*vec[maxi]/(nc);
                        }
                    }
                    if(ui->viewBox->currentIndex()==0){ // median
                        for (size_t p = 0; p != idx.size(); ++p) {idx[p] = p; }
                        sort (idx.begin (), idx.end (), compare_index<vector<unsigned char> &>(vec));

                        mini = idx[0];
                        maxi = idx[nimg-1];
                        if((vec.size()/2) % 2 == 0)
                            medi = idx[(vec.size())/2 + 1];
                        else
                            medi = idx[(vec.size())/2];

                        for(int k=0;k<3;k++){
                            Vec3b col=image.at<Vec3b>(j,i);
                            float nc = col[0]+col[1]+col[2];
                            val[k] = 3*col[k]*vec[medi]/(nc);
                        }

                    }
                    else if(ui->viewBox->currentIndex()==1){ //average
                        for (size_t p = 0; p != idx.size(); ++p) {avg=avg+vec[p];}
                        avg /= vec.size();

                        for(int k=0;k<3;k++){
                            Vec3b col=image.at<Vec3b>(j,i);
                            float nc = col[0]+col[1]+col[2];
                            val[k] = 3*col[k]*avg/(nc);
                        }

                    }
                    else if(ui->viewBox->currentIndex()==2){
                        // interpolate

                        Vec3b col=image.at<Vec3b>(j,i);
                        float nc = col[0]+col[1]+col[2];

                        int nen=7;
                        if (flagInterp == false) {
                            rbfcreate(2, 1, modelInterpolation);
                            arrayInterpolation.setlength(nen, 3);

                            /*
                    for(int k=0;k<vec.size();k++){
                        float dv=(lx-dirs[k][0])*(lx-dirs[k][0])+(ly-dirs[k][1])*(ly-dirs[k][1])+(lz-dirs[k][2])*(lz-dirs[k][2]);
                        dist.push_back(dv);
                    }

                    iv.resize(vec.size());
                    for (size_t p = 0; p != iv.size(); ++p) iv[p] = p;
                    sort (iv.begin (), iv.end (), compare_index<vector<float> &>(dist));
                    dist.clear();*/

                            for(int k1=0; k1<nen; k1++)
                                for(int k2=0; k2<2; k2++) {
                                    arrayInterpolation(k1, k2) = dirs[iv[k1]][k2];
                                }

                            flagInterp=true;
                        }

                        for(int k=0; k<nen; k++)
                            arrayInterpolation(k, 2) = vec[iv[k]];

                        rbfsetpoints(modelInterpolation, arrayInterpolation);
                        float rr = ui->rbfSpinBox->value();
                        rbfsetalgomultilayer(modelInterpolation, rr, 1, 10e-3);

                        rbfbuildmodel(modelInterpolation, repInterpolation);

                        double res = rbfcalc2(modelInterpolation, lx, ly);

                        for (int k=0; k<3; k++)
                            val[k] = std::max(0.0,std::min(255.0,3*col[k]*res/(nc)));


                        /*

                    float distvector[] = {(lx-dirs[iv[0]][0])*(lx-dirs[iv[0]][0])+(ly-dirs[iv[0]][1])*(ly-dirs[iv[0]][1])+(lz-dirs[iv[0]][2])*(lz-dirs[iv[0]][2]),
                                          (lx-dirs[iv[1]][0])*(lx-dirs[iv[1]][0])+(ly-dirs[iv[1]][1])*(ly-dirs[iv[1]][1])+(lz-dirs[iv[1]][2])*(lz-dirs[iv[1]][2]),
                                          (lx-dirs[iv[2]][0])*(lx-dirs[iv[2]][0])+(ly-dirs[iv[2]][1])*(ly-dirs[iv[2]][1])+(lz-dirs[iv[2]][2])*(lz-dirs[iv[2]][2])};
                    float distf = distvector[0]+distvector[1]+distvector[2];

                    val[0] = 0;
                    val[1] = 0;
                    val[2] = 0;

                    for(int k=0;k<3;k++){
                        Vec3b col=image.at<Vec3b>(j,i);
                        float nc = col[0]+col[1]+col[2];
                       // float pv =(vec[iv[0]]*dist[iv[0]] + vec[iv[1]]*dist[iv[0]])/(dist[iv[0]]+dist[iv[1]]);
                       // vec[iv[0]] è il nearest neighbor il resto sarebbe per
                       // convertire a 8 bit e applicare la tinta (che forse non è corretta, ma non te ne preoccupare)
                        for (int i = 0; i < 3; i++)
                            val[k] = val[k] + (3*col[k]*vec[iv[i]]/(nc*256))*((distf-distvector[i])/distf);
                    }*/
                    }
                    else if(ui->viewBox->currentIndex()==4){
                        // interpolate removing an image on the missing image directions

                        Vec3b col=image.at<Vec3b>(j,i);
                        float nc = col[0]+col[1]+col[2];

                        int nen=7;
                        if (flagInterp == false) {
                            rbfcreate(2, 1, modelInterpolation);
                            arrayInterpolation.setlength(nen, 3);

                            qDebug() << dirs[0][0]  << dirs[0][1]  << dirs[0][2];
                            qDebug() << dirs[1][0]  << dirs[1][1]  << dirs[1][2];



                            // prendo i vicini escluso il primo
                            for(int k1=0; k1<nen; k1++)
                                for(int k2=0; k2<2; k2++) {
                                    arrayInterpolation(k1, k2) = dirs[iv[k1+1]][k2];
                                }

                            flagInterp=true;
                        }
                        lx=dirs[iv[0]][0];
                        ly=dirs[iv[0]][1];

                        for(int k=0; k<nen; k++)
                            arrayInterpolation(k, 2) = vec[iv[k+1]];

                        rbfsetpoints(modelInterpolation, arrayInterpolation);
                        float rr = ui->rbfSpinBox->value();
                        rbfsetalgomultilayer(modelInterpolation, rr, 1, 10e-3);

                        rbfbuildmodel(modelInterpolation, repInterpolation);

                        double res = rbfcalc2(modelInterpolation, lx, ly);

                        for (int k=0; k<3; k++)
                            val[k] = 3*col[k]*res/(nc);
                        for (int k=0; k<3; k++)
                            valcomp[k] = 3*col[k]*vec[iv[0]]/nc;

                        comparedImg.at<Vec3b>(j,i) = valcomp;
                    }


                    shownImg.at<Vec3b>(j,i) = val;


                }
                // 16 bit LRGB
                if(type==2){

                    // read pixel data and do some processing. on this skeleton we can
                    // develop fitters, estimate normals, detect shadows and edges, etc.

                    in.readRawData((char*)&vals[0], nimg*2);
                    vector<unsigned short> vec (&vals[0], &vals[0]+nimg );
                    int medi, mini, maxi;

                    int out;

                    float avg;
                    /*  vector<size_t> idx(vec.size());
                for (size_t p = 0; p != idx.size(); ++p) {idx[p] = p; avg=avg+vec[p];}
                sort (idx.begin (), idx.end (), compare_index<vector<unsigned short> &>(vec));
                avg /= vec.size();

                mini = idx[0];
                maxi = idx[nimg-1];
                if((vec.size()/2) % 2 == 0)
                    medi = idx[(vec.size())/2 + 1];
                else
                    medi = idx[(vec.size())/2];
*/

                    /*  Vec3f cc; cc[0] = 0.2126; cc[1] = 0.7152; cc[2] = 0.0722;*/

                    if(ui->viewBox->currentIndex()==0){
                        vector<size_t> idx(vec.size());
                        for (size_t p = 0; p != idx.size(); ++p) {idx[p] = p;}
                        sort (idx.begin (), idx.end (), compare_index<vector<unsigned short> &>(vec));

                        mini = idx[0];
                        maxi = idx[nimg-1];
                        if((vec.size()/2) % 2 == 0)
                            medi = idx[(vec.size())/2 + 1];
                        else
                            medi = idx[(vec.size())/2];

                        for(int k=0;k<3;k++){
                            Vec3b col=image.at<Vec3b>(j,i);

                            float nc = col[0]+col[1]+col[2];
                            val[k] = 3*col[k]*vec[medi]/(nc*256);
                        }
                    }
                    else if(ui->viewBox->currentIndex()==1){
                        vector<size_t> idx(vec.size());
                        for (size_t p = 0; p != idx.size(); ++p) {avg=avg+vec[p];}
                        avg /= vec.size();

                        for(int k=0;k<3;k++){
                            Vec3b col=image.at<Vec3b>(j,i);

                            float nc = col[0]+col[1]+col[2];

                            val[k] = 3*col[k]*avg/(nc*255);

                        }
                    }
                    else if(ui->viewBox->currentIndex()==2){
                        // interpolation

                        Vec3b col=image.at<Vec3b>(j,i);
                        float nc = col[0]+col[1]+col[2];
                        int nen=7;
                        if (flagInterp == false) {
                            rbfcreate(2, 1, modelInterpolation);
                            arrayInterpolation.setlength(nen, 3);

                            for(int k=0;k<vec.size();k++){
                                float dv=(lx-dirs[k][0])*(lx-dirs[k][0])+(ly-dirs[k][1])*(ly-dirs[k][1])+(lz-dirs[k][2])*(lz-dirs[k][2]);
                                dist.push_back(dv);
                            }


                            for(int k1=0; k1<nen; k1++)
                                for(int k2=0; k2<2; k2++) {
                                    arrayInterpolation(k1, k2) = dirs[iv[k1]][k2];
                                }

                            flagInterp=true;
                        }

                        for(int k=0; k<nen; k++)
                            arrayInterpolation(k, 2) = vec[iv[k]];


                        rbfsetpoints(modelInterpolation, arrayInterpolation);
                        float rr = ui->rbfSpinBox->value();
                        rbfsetalgomultilayer(modelInterpolation,rr, 1, 10e-3);

                        rbfbuildmodel(modelInterpolation, repInterpolation);

                        double res = rbfcalc2(modelInterpolation, lx, ly);

                        for (int k=0; k<3; k++)
                            val[k] = std::min(255.0,std::max(0.0,3*col[k]*res/(nc*255)));

                        /*//ORA BANALMENTE PRENDE LA DISTANZA, ORDINA e prende nearst neighbor
                    for(int k=0;k<vec.size();k++){
                        float dv=(lx-dirs[k][0])*(lx-dirs[k][0])+(ly-dirs[k][1])*(ly-dirs[k][1])+(lz-dirs[k][2])*(lz-dirs[k][2]);
                        dist.push_back(dv);
                    }

                    vector<size_t> iv(vec.size());
                    for (size_t p = 0; p != iv.size(); ++p) iv[p] = p;
                    sort (iv.begin (), iv.end (), compare_index<vector<float> &>(dist));
                    dist.clear();

                    float distvector[] = {(lx-dirs[iv[0]][0])*(lx-dirs[iv[0]][0])+(ly-dirs[iv[0]][1])*(ly-dirs[iv[0]][1])+(lz-dirs[iv[0]][2])*(lz-dirs[iv[0]][2]),
                                          (lx-dirs[iv[1]][0])*(lx-dirs[iv[1]][0])+(ly-dirs[iv[1]][1])*(ly-dirs[iv[1]][1])+(lz-dirs[iv[1]][2])*(lz-dirs[iv[1]][2]),
                                          (lx-dirs[iv[2]][0])*(lx-dirs[iv[2]][0])+(ly-dirs[iv[2]][1])*(ly-dirs[iv[2]][1])+(lz-dirs[iv[2]][2])*(lz-dirs[iv[2]][2])};
                    float distf = distvector[0]+distvector[1]+distvector[2];

                    val[0] = 0;
                    val[1] = 0;
                    val[2] = 0;

                    for(int k=0;k<3;k++){
                        Vec3b col=image.at<Vec3b>(j,i);
                        float nc = col[0]+col[1]+col[2];
                       // float pv =(vec[iv[0]]*dist[iv[0]] + vec[iv[1]]*dist[iv[0]])/(dist[iv[0]]+dist[iv[1]]);
                       // vec[iv[0]] è il nearest neighbor il resto sarebbe per
                       // convertire a 8 bit e applicare la tinta (che forse non è corretta, ma non te ne preoccupare)
                        for (int i = 0; i < 3; i++)
                            val[k] = val[k] + (3*col[k]*vec[iv[i]]/(nc*256))*((distf-distvector[i])/distf);
                    }*/
                    }

                    shownImg.at<Vec3b>(j,i) = val;
                }

            }
    else if(type>2){
        // 8 bit RGB
        if(type==3)
            for(int cc=0; cc<3;cc++) { // loop over colors
                for(int j=0;j<size[1];j++)
                    for(int i=0;i<size[0];i++)
                    {
                        if(dirtype==2) // interpolated dirs: estimate pixel-specific direction
                            for(int k=0;k<nimg;k++)
                            {
                                dirs[k][0]=dircoeffs[k][0]*i+dircoeffs[k][1]*j+dircoeffs[k][2];
                                dirs[k][1]=dircoeffs[k][3]*i+dircoeffs[k][4]*j+dircoeffs[k][5];
                                dirs[k][2]=dircoeffs[k][6]*i+dircoeffs[k][7]*j+dircoeffs[k][8];

                                if(i==0 && j==0)
                                    qDebug() << " 0,0 - " << dirs[k][0] << " " << dirs[k][1] << " "  << dirs[k][2];
                                if(i==size[0]-1 && j==size[1]-1)
                                    qDebug() << " end  - " << dirs[k][0] << " " << dirs[k][1] << " "  << dirs[k][2];
                            }



                        in.readRawData((char*)&valc[0], nimg);
                        vector<unsigned char> vec (&valc[0], &valc[0]+nimg );
                        int medi, mini, maxi;
                        vector<int> shad;
                        vector<int> shad2, tmpv;
                        // 8 bit version,

                        // read pixel data and do some processing. on this skeleton we can
                        // develop fitters, estimate normals, detect shadows and edges, etc.

                        float avg;
                        vector<size_t> idx(vec.size());
                        /* for (size_t p = 0; p != idx.size(); ++p) {idx[p] = p; avg=avg+vec[p];}
            sort (idx.begin (), idx.end (), compare_index<vector<unsigned char> &>(vec));
            avg /= vec.size();

            mini = idx[0];
            maxi = idx[nimg-1];
            if((vec.size()/2) % 2 == 0)
                medi = idx[(vec.size())/2 + 1];
            else
                medi = idx[(vec.size())/2];

*/
                        /*  Vec3f cc; cc[0] = 0.2126; cc[1] = 0.7152; cc[2] = 0.0722;*/

                        if(ui->viewBox->currentIndex()==3){ // max
                            for (size_t p = 0; p != idx.size(); ++p) {idx[p] = p; }
                            sort (idx.begin (), idx.end (), compare_index<vector<unsigned char> &>(vec));
                            mini = idx[0];
                            maxi = idx[nimg-1];
                            shownImg.at<Vec3b>(j,i)[cc] = vec[maxi];
                        }
                        if(ui->viewBox->currentIndex()==0){ // median

                            for (size_t p = 0; p != idx.size(); ++p) {idx[p] = p; }
                            sort (idx.begin (), idx.end (), compare_index<vector<unsigned char> &>(vec));

                            mini = idx[0];
                            maxi = idx[nimg-1];
                            if((vec.size()/2) % 2 == 0)
                                medi = idx[(vec.size())/2 + 1];
                            else
                                medi = idx[(vec.size())/2];

                            shownImg.at<Vec3b>(j,i)[cc] = vec[medi];


                        }
                        else if(ui->viewBox->currentIndex()==1){ // mean
                            for (size_t p = 0; p != idx.size(); ++p) { avg=avg+vec[p];}
                            avg /= vec.size();

                            shownImg.at<Vec3b>(j,i)[cc] = avg;

                        }
                        else if(ui->viewBox->currentIndex()==2){ // directional
                            // interpolate
                            // QUI DA IMPLEMENTARE!

                            //  qDebug() << i << " " << j  << "\n";

                            /*              float lx = ui->lxSpinBox->value();
                float ly = ui->lySpinBox->value();
                float lz = (1-lx*lx-ly*ly);*/

                            // Vec3b col=image.at<Vec3b>(j,i);
                            // float nc = col[0]+col[1]+col[2];
                            int nen=7;

                            if (flagInterp == false) {
                                rbfcreate(2, 1, modelInterpolation);
                                arrayInterpolation.setlength(nen, 3);

                                /*
                for(int k=0;k<vec.size();k++){
                    float dv=(lx-dirs[k][0])*(lx-dirs[k][0])+(ly-dirs[k][1])*(ly-dirs[k][1])+(lz-dirs[k][2])*(lz-dirs[k][2]);
                    dist.push_back(dv);
                }

               iv.resize(vec.size());
                for (size_t p = 0; p != iv.size(); ++p) iv[p] = p;
                sort (iv.begin (), iv.end (), compare_index<vector<float> &>(dist));
                dist.clear();
*/

                                for(int k1=0; k1<nen; k1++)
                                    for(int k2=0; k2<2; k2++) {
                                        arrayInterpolation(k1, k2) = dirs[iv[k1]][k2];
                                    }

                                flagInterp=true;
                            }

                            for(int k=0; k<nen; k++)
                                arrayInterpolation(k, 2) = vec[iv[k]];


                            rbfsetpoints(modelInterpolation, arrayInterpolation);
                            float rr = ui->rbfSpinBox->value();
                            rbfsetalgomultilayer(modelInterpolation, rr, 1, 10e-3);

                            rbfbuildmodel(modelInterpolation, repInterpolation);

                            double res = rbfcalc2(modelInterpolation, lx, ly);

                            // for (int k=0; k<3; k++)
                            //     val[cc] = 3*col[cc]*res/(nc);

                            shownImg.at<Vec3b>(j,i)[cc] = std::max(0.0,std::min(255.0,res));

                            /*//ORA BANALMENTE PRENDE LA DISTANZA, ORDINA e prende nearst neighbor
                for(int k=0;k<vec.size();k++){
                    float dv=(lx-dirs[k][0])*(lx-dirs[k][0])+(ly-dirs[k][1])*(ly-dirs[k][1])+(lz-dirs[k][2])*(lz-dirs[k][2]);
                    dist.push_back(dv);
                }

                vector<size_t> iv(vec.size());
                for (size_t p = 0; p != iv.size(); ++p) iv[p] = p;
                sort (iv.begin (), iv.end (), compare_index<vector<float> &>(dist));
                dist.clear();

                float distvector[] = {(lx-dirs[iv[0]][0])*(lx-dirs[iv[0]][0])+(ly-dirs[iv[0]][1])*(ly-dirs[iv[0]][1])+(lz-dirs[iv[0]][2])*(lz-dirs[iv[0]][2]),
                                      (lx-dirs[iv[1]][0])*(lx-dirs[iv[1]][0])+(ly-dirs[iv[1]][1])*(ly-dirs[iv[1]][1])+(lz-dirs[iv[1]][2])*(lz-dirs[iv[1]][2]),
                                      (lx-dirs[iv[2]][0])*(lx-dirs[iv[2]][0])+(ly-dirs[iv[2]][1])*(ly-dirs[iv[2]][1])+(lz-dirs[iv[2]][2])*(lz-dirs[iv[2]][2])};
                float distf = distvector[0]+distvector[1]+distvector[2];

                val[0] = 0;
                val[1] = 0;
                val[2] = 0;

                for(int k=0;k<3;k++){
                    Vec3b col=image.at<Vec3b>(j,i);
                    float nc = col[0]+col[1]+col[2];
                   // float pv =(vec[iv[0]]*dist[iv[0]] + vec[iv[1]]*dist[iv[0]])/(dist[iv[0]]+dist[iv[1]]);
                   // vec[iv[0]] è il nearest neighbor il resto sarebbe per
                   // convertire a 8 bit e applicare la tinta (che forse non è corretta, ma non te ne preoccupare)
                    for (int i = 0; i < 3; i++)
                        val[k] = val[k] + (3*col[k]*vec[iv[i]]/(nc*256))*((distf-distvector[i])/distf);
                }*/
                        }
                        else if(ui->viewBox->currentIndex()==4){
                            // interpolate removing an image on the missing image directions

                            int nen=7;
                            if (flagInterp == false) {
                                rbfcreate(2, 1, modelInterpolation);
                                arrayInterpolation.setlength(nen, 3);

                                qDebug() << dirs[0][0]  << dirs[0][1]  << dirs[0][2];
                                qDebug() << dirs[1][0]  << dirs[1][1]  << dirs[1][2];

                                // prendo i vicini escluso il primo
                                for(int k1=0; k1<nen; k1++)
                                    for(int k2=0; k2<2; k2++) {
                                        arrayInterpolation(k1, k2) = dirs[iv[k1+1]][k2];
                                    }

                                flagInterp=true;
                            }
                            lx=dirs[iv[0]][0];
                            ly=dirs[iv[0]][1];

                            for(int k=0; k<nen; k++)
                                arrayInterpolation(k, 2) = vec[iv[k+1]];

                            rbfsetpoints(modelInterpolation, arrayInterpolation);
                            float rr = ui->rbfSpinBox->value();
                            rbfsetalgomultilayer(modelInterpolation, rr, 1, 10e-3);

                            rbfbuildmodel(modelInterpolation, repInterpolation);

                            double res = rbfcalc2(modelInterpolation, lx, ly);

                            shownImg.at<Vec3b>(j,i)[cc] = std::min(255.0,std::max(0.0,res));


                            comparedImg.at<Vec3b>(j,i)[cc] = vec[iv[0]];
                        }


                    }

            } // end loop over colors

        if(type==4) // 16 bit RGB
            for(int cc=0; cc<3;cc++) { // loop over colors
                for(int j=0;j<size[1];j++)
                    for(int i=0;i<size[0];i++)
                    {
                        if(dirtype==2) // interpolated dirs: estimate pixel-specific direction
                            for(int k=0;k<nimg;k++)
                            {
                                dirs[k][0]=dircoeffs[k][0]*i+dircoeffs[k][1]*j+dircoeffs[k][2];
                                dirs[k][1]=dircoeffs[k][3]*i+dircoeffs[k][4]*j+dircoeffs[k][5];
                                dirs[k][2]=dircoeffs[k][6]*i+dircoeffs[k][7]*j+dircoeffs[k][8];

                                if(i==0 && j==0)
                                    qDebug() << " 0,0 - " << dirs[k][0] << " " << dirs[k][1] << " "  << dirs[k][2];
                                if(i==size[0]-1 && j==size[1]-1)
                                    qDebug() << " end  - " << dirs[k][0] << " " << dirs[k][1] << " "  << dirs[k][2];
                            }

                        in.readRawData((char*)&vals[0], nimg*2);
                        vector<unsigned short> vec (&vals[0], &vals[0]+nimg );
                        int medi, mini, maxi;
                        vector<int> shad;
                        vector<int> shad2, tmpv;
                        // 16 bit version,

                        // read pixel data and do some processing. on this skeleton we can
                        // develop fitters, estimate normals, detect shadows and edges, etc.

                        float avg;
                        vector<size_t> idx(vec.size());


                        if(ui->viewBox->currentIndex()==3){ // max
                            for (size_t p = 0; p != idx.size(); ++p) {idx[p] = p; }
                            sort (idx.begin (), idx.end (), compare_index<vector<unsigned short> &>(vec));
                            mini = idx[0];
                            maxi = idx[nimg-1];
                            shownImg.at<Vec3b>(j,i)[cc] = vec[maxi]/256.0;
                        }
                        if(ui->viewBox->currentIndex()==0){ // median

                            for (size_t p = 0; p != idx.size(); ++p) {idx[p] = p; }
                            sort (idx.begin (), idx.end (), compare_index<vector<unsigned short> &>(vec));

                            mini = idx[0];
                            maxi = idx[nimg-1];
                            if((vec.size()/2) % 2 == 0)
                                medi = idx[(vec.size())/2 + 1];
                            else
                                medi = idx[(vec.size())/2];

                            shownImg.at<Vec3b>(j,i)[cc] = vec[medi]/256.0;


                        }
                        else if(ui->viewBox->currentIndex()==1){ // mean
                            for (size_t p = 0; p != idx.size(); ++p) { avg=avg+vec[p];}
                            avg /= vec.size();

                            shownImg.at<Vec3b>(j,i)[cc] = avg/256.0;

                        }
                        else if(ui->viewBox->currentIndex()==2){ // directional

                            int nen=7;

                            if (flagInterp == false) {
                                rbfcreate(2, 1, modelInterpolation);
                                arrayInterpolation.setlength(nen, 3);


                                for(int k1=0; k1<nen; k1++)
                                    for(int k2=0; k2<2; k2++) {
                                        arrayInterpolation(k1, k2) = dirs[iv[k1]][k2];
                                    }

                                flagInterp=true;
                            }

                            for(int k=0; k<nen; k++)
                                arrayInterpolation(k, 2) = vec[iv[k]];


                            rbfsetpoints(modelInterpolation, arrayInterpolation);
                            float rr = ui->rbfSpinBox->value();
                            rbfsetalgomultilayer(modelInterpolation, rr, 1, 10e-3);

                            rbfbuildmodel(modelInterpolation, repInterpolation);

                            double res = rbfcalc2(modelInterpolation, lx, ly);


                            shownImg.at<Vec3b>(j,i)[cc] = std::min(255.0,res/256.0);

                            /*//ORA BANALMENTE PRENDE LA DISTANZA, ORDINA e prende nearst neighbor
                for(int k=0;k<vec.size();k++){
                    float dv=(lx-dirs[k][0])*(lx-dirs[k][0])+(ly-dirs[k][1])*(ly-dirs[k][1])+(lz-dirs[k][2])*(lz-dirs[k][2]);
                    dist.push_back(dv);
                }

                vector<size_t> iv(vec.size());
                for (size_t p = 0; p != iv.size(); ++p) iv[p] = p;
                sort (iv.begin (), iv.end (), compare_index<vector<float> &>(dist));
                dist.clear();

                float distvector[] = {(lx-dirs[iv[0]][0])*(lx-dirs[iv[0]][0])+(ly-dirs[iv[0]][1])*(ly-dirs[iv[0]][1])+(lz-dirs[iv[0]][2])*(lz-dirs[iv[0]][2]),
                                      (lx-dirs[iv[1]][0])*(lx-dirs[iv[1]][0])+(ly-dirs[iv[1]][1])*(ly-dirs[iv[1]][1])+(lz-dirs[iv[1]][2])*(lz-dirs[iv[1]][2]),
                                      (lx-dirs[iv[2]][0])*(lx-dirs[iv[2]][0])+(ly-dirs[iv[2]][1])*(ly-dirs[iv[2]][1])+(lz-dirs[iv[2]][2])*(lz-dirs[iv[2]][2])};
                float distf = distvector[0]+distvector[1]+distvector[2];

                val[0] = 0;
                val[1] = 0;
                val[2] = 0;

                for(int k=0;k<3;k++){
                    Vec3b col=image.at<Vec3b>(j,i);
                    float nc = col[0]+col[1]+col[2];
                   // float pv =(vec[iv[0]]*dist[iv[0]] + vec[iv[1]]*dist[iv[0]])/(dist[iv[0]]+dist[iv[1]]);
                   // vec[iv[0]] è il nearest neighbor il resto sarebbe per
                   // convertire a 8 bit e applicare la tinta (che forse non è corretta, ma non te ne preoccupare)
                    for (int i = 0; i < 3; i++)
                        val[k] = val[k] + (3*col[k]*vec[iv[i]]/(nc*256))*((distf-distvector[i])/distf);
                }*/
                        }
                    }

            } // end loop over colors

    }


    if(ui->viewBox->currentIndex()==0){
        ui->msgBox->setText("Median image");
    }
    if(ui->viewBox->currentIndex()==1){
        ui->msgBox->setText("Median image");
    }
    if(ui->viewBox->currentIndex()==2){
        ui->msgBox->setText("Rbf interpolated image");
        ui->msgBox->append("Relighted from" + QString::number(lx) +  " " + QString::number(ly));
    }
    if(ui->viewBox->currentIndex()==3){
        ui->msgBox->setText("max image");
    }
    iw->imageLabel->setPixmap(QPixmap::fromImage(QImage(shownImg.data,shownImg.cols,shownImg.rows,shownImg.step,QImage::Format_RGB888)));

    filed.close();
    cv::cvtColor(shownImg,shownImg, cv::COLOR_BGR2RGB);


    if(ui->viewBox->currentIndex()==4)
    { cv::cvtColor(comparedImg,comparedImg, cv::COLOR_BGR2RGB);
       QString name1 = "original_" + QString::number(ui->spinBox->value()) + ".png";
        imwrite(name1.toStdString(),comparedImg);
      name1 = "relightedRBF_"  + QString::number(ui->rbfSpinBox->value()) + "_" + QString::number(ui->spinBox->value())  + ".png";

        imwrite(name1.toStdString(),shownImg);
    }
    else
     imwrite("result.png",shownImg);
}

void apTool::on_pushButton_6_clicked()
{
    iw->white1->setGeometry(QRect(0,0,0,0));
    iw->active = 11;
}

void apTool::on_pushButton_7_clicked()
{
    iw->white2->setGeometry(QRect(0,0,0,0));
    iw->active = 12;
}

void apTool::on_pushButton_2_clicked()
{
    iw->CoordinatesSet1=false;
    iw->point1->setGeometry(QRect(0,0,0,0));
    iw->active=1;
}

void apTool::on_pushButton_3_clicked()
{
    iw->CoordinatesSet2=false;
    iw->point2->setGeometry(QRect(0,0,0,0));
    iw->active=2;
}

void apTool::on_clearButton_clicked()
{
    iw->white1->setGeometry(QRect(0,0,0,0));
    iw->point1->setGeometry(QRect(0,0,0,0));
    iw->CoordinatesSet1=false;
}

void apTool::on_clearButton_2_clicked()
{
    iw->white2->setGeometry(QRect(0,0,0,0));
    iw->point2->setGeometry(QRect(0,0,0,0));
    iw->CoordinatesSet2=false;
}

void apTool::on_lWidget_clicked(int lx, int ly) {

    double flx = 0.00;
    double fly = 0.00;

    if (lx < (ui->lWidget->width()/(double)2))
        flx = -((ui->lWidget->width()/(double)2)-lx)/(double)((ui->lWidget->width()/(double)2));
    else if (lx > (ui->lWidget->width()/(double)2))
        flx = (lx-(ui->lWidget->width()/(double)2))/(double)((ui->lWidget->width()/(double)2));

    if (ly < (ui->lWidget->height()/(double)2))
        fly = ((ui->lWidget->height()/(double)2)-ly)/(double)(ui->lWidget->height()/(double)2);
    else if (ly > (ui->lWidget->height()/(double)2))
        fly = -(ly-(ui->lWidget->height()/(double)2))/(double)(ui->lWidget->height()/(double)2);

    ui->lxSpinBox->setValue(flx);
    ui->lySpinBox->setValue(fly);

    QPixmap image = drawLWidget();
    QPainter painter(&image);

    painter.setBrush(QBrush(Qt::yellow));
    QRectF rectangle(lx-(5), ly-(5), 10.0, 10.0);
    painter.drawEllipse(rectangle);
    ui->lWidget->setPixmap(image);
}

void apTool::on_pushButton_4_clicked()
{
    showAP("AP 1");
}

void apTool::on_pushButton_5_clicked()
{
    showAP("AP 2");
}

void apTool::showAP(String windowName)
{

    QString filename = ui->fileNameLine->text();
    QFile file(filename);
    int last= filename.lastIndexOf(QDir::separator());
    QString folder=filename.left(last+1);

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    double zoom = iw->scaleFactor;
    int x;
    int y;
    int wid=0;
    int hei=0;

    if (windowName.compare("AP 1") == 0) {
        if (iw->CoordinatesSet1==true) {
            x = iw->x1/zoom;
            y = iw->y1/zoom;
        }
        else {
            x = iw->white1->x()/zoom;
            y = iw->white1->y()/zoom;
            wid = iw->white1->width();
            hei = iw->white1->height();
        }
    }
    else if (windowName.compare("AP 2") == 0) {
        if (iw->CoordinatesSet2==true) {
            x = iw->x2/zoom;
            y = iw->y2/zoom;
        }
        else {
            x = iw->white2->x()/zoom;
            y = iw->white2->y()/zoom;
            wid = iw->white2->width();
            hei = iw->white2->height();
        }
    }
    else return;

    if (iw->CoordinatesSet1==false && iw->CoordinatesSet2==false && (wid == 0 || hei ==0))
        return;

    int type=0;
    int dirtype=0;
    int size[2]={0,0};
    int nimg=0;
    Vec3f* dirs;
    float** dircoeffs;
    QString chroma_img;


    // Reading header file

    QTextStream textStream(&file);

    QString line = textStream.readLine();

    if (line.isNull())
        return;

    QStringList parts = line.split(" ");
    if(parts[0] != QString("APA"))
        return;

    while(!(line.isNull())){
        line = textStream.readLine();
        // qDebug() << line;

        parts = line.split(" ");
        if(parts[0]== "LUMINANCE_TYPE" && parts[1]== "UNSIGNED_CHAR" ){
            qDebug() << "UNSIGNED CHAR";
            type=1;
        }
        if(parts[0]== "LUMINANCE_TYPE" && parts[1]== "UNSIGNED_SHORT" ){
            qDebug() << "UNSIGNED SHORT";
            type=2;
        }

        if(parts[0]== "COLOR_TYPE" && parts[1]== "UNSIGNED_CHAR" ){
            qDebug() << "COLOR UNSIGNED CHAR";
            type=3;
        }
        if(parts[0]== "COLOR_TYPE" && parts[1]== "UNSIGNED_SHORT" ){
            qDebug() << "COLOR UNSIGNED SHORT";
            type=4;
        }



        if(parts[0]== "CHROMA_IMAGE"  ){
            //chroma_img = parts[1];
            chroma_img = folder+parts[1];

            qDebug() << chroma_img;
        }

        if(parts[0]== "IMAGE_SIZE"  ){
            size[0]=parts[1].toInt();
            size[1]=parts[2].toInt();
            qDebug() << "size " << size[0] << " " << size[1];
        }

        if(parts[0]== "N_IMAGES"){
            nimg = parts[1].toInt();
            qDebug() << "number of images " <<  nimg;

        }
        if(parts[0]== "DIR_TYPE"){
            if(parts[1] == "constant"){
                dirtype=1;
                qDebug() << "Dir type constant";}
            if(parts[1] == "interpolated"){
                dirtype=2;
                qDebug() << "Dir type interpolated";}
        }


        if(parts[0]== "DIRECTION_COEFFICIENTS"){
            if(dirtype==2){
                dirs = new Vec3f[nimg];
                dircoeffs = new float*[nimg];
                for(int i=0;i<nimg;i++)   {
                    line = textStream.readLine();
                    dircoeffs[i]=new float[9];

                    if (line.isNull()){
                        qDebug() << "error";
                        return;
                    }
                    parts = line.split(" ");
                    if(parts.size() < 6){
                        qDebug() << "error";
                        return;
                    }

                    for(int j=0;j<9;j++)
                        dircoeffs[i][j] = parts[j].toFloat();

                    //qDebug() << dircoeffs[i][0] << " " << dircoeffs[i][1] << " " << dircoeffs[i][5];
                    float c_x=size[0]/2;float c_y=size[1]/2;
                    dirs[i][0]=dircoeffs[i][0]*c_x+dircoeffs[i][1]*c_y+dircoeffs[i][2];
                    dirs[i][1]=dircoeffs[i][3]*c_x+dircoeffs[i][4]*c_y+dircoeffs[i][5];
                    dirs[i][2]=dircoeffs[i][6]*c_x+dircoeffs[i][7]*c_y+dircoeffs[i][8];
                    qDebug() << dirs[i][0] << " " << dirs[i][1] << " "  << dirs[i][2];

                    // computing interpolated direction in the center
                    //
                }
            }
        }


        //if light directions are constant
        if(parts[0] == "LIGHT_DIRECTIONS"){
            if(dirtype==1){
                dirs = new Vec3f[nimg];
                for(int i=0;i<nimg;i++)   {
                    line = textStream.readLine();
                    if (line.isNull()){
                        qDebug() << "error";
                        return;
                    }
                    parts = line.split(" ");
                    if(parts.size() < 3){
                        qDebug() << "error";
                        return;
                    }

                    dirs[i][0]=parts[0].toFloat();
                    dirs[i][1]=parts[1].toFloat();
                    dirs[i][2]=parts[2].toFloat();

                }
            }

        }
    }

    file.close();

    cv::Mat image;
    image = cv::imread(chroma_img.toStdString(), CV_LOAD_IMAGE_COLOR);

    filename.replace(".aph",".apd");
    QFile filed(filename);

    if (!filed.open(QIODevice::ReadOnly))
        return;

    unsigned short* vals = new unsigned short[nimg];
    unsigned char* valc = new unsigned char[nimg];

    QDataStream in(&filed);

    if(type==1){
        for(int j=1;j<y;j++)
            for (int i=1;i<=size[0];i++)
                in.readRawData((char*)&valc[0], nimg);

        for(int i=1;i<=x;i++)
            in.readRawData((char*)&valc[0], nimg);
    }

    if(type==2){
        for(int j=1;j<y;j++)
            for (int i=1;i<=size[0];i++)
                in.readRawData((char*)&vals[0], nimg*2);

        for(int i=1;i<=x;i++)
            in.readRawData((char*)&vals[0], nimg*2);
    }


    int scale = 2;
    Mat shownImg(scale*101,scale*101,CV_8UC3);

    rbfmodel modelInterpolation;
    rbfreport repInterpolation;
    real_2d_array arrayInterpolation;
    Vec3f val;
    Vec3b col;
    float nc = 0;

    rbfcreate(2, 1, modelInterpolation);
    arrayInterpolation.setlength(nimg, 3);


    vector<unsigned char> vecc;
    vector<unsigned short> vecs;


    for(int k1=0; k1<nimg; k1++)
        for(int k2=0; k2<2; k2++) {
            arrayInterpolation(k1, k2) = dirs[k1][k2];
        }

    if (wid == 0 || hei ==0){

        if(type==2){
            vector<unsigned short> vecs (&vals[0], &vals[0]+nimg);

            col=image.at<Vec3b>(y,x);
            nc = col[0] + col[1] + col[2];

            for(int k=0; k<nimg; k++)
                arrayInterpolation(k, 2) = vecs[k];
        }

        if(type==1){
            vector<unsigned char> vecc (&valc[0], &valc[0]+nimg);

            col=image.at<Vec3b>(y,x);
            nc = col[0] + col[1] + col[2];

            for(int k=0; k<nimg; k++)
                arrayInterpolation(k, 2) = vecc[k];
        }
    }

    else {

        float conv1[3] = {0,0,0};
        float conv2[nimg];
        for(int k=0; k<nimg; k++)
            conv2[k] = 0;
        for(int k=0; k<3; k++)
            col[k] = 0;


        if(type==2)
            for (int j = 0; j < hei; j++) {
                for (int i = 0; i < wid; i++) {

                    Vec3b colPixel=image.at<Vec3b>(j,i);
                    for(int k=0; k<3; k++)
                        conv1[k] += colPixel[k];
                    nc += colPixel[0]+colPixel[1]+colPixel[2];
                    in.readRawData((char*)&vals[0], nimg*2);
                    vecs.assign(&vals[0], &vals[0]+nimg);
                    for(int k=0; k<nimg; k++)
                        conv2[k] += vecs[k];

                }
                for(int k=0; k<=size[0]-wid; k++)
                    in.readRawData((char*)&vals[0], nimg*2);
                vecs.assign(&vals[0], &vals[0]+nimg);
            }
        if(type==1)
            for (int j = 0; j < hei; j++) {
                for (int i = 0; i < wid; i++) {

                    Vec3b colPixel=image.at<Vec3b>(j,i);
                    for(int k=0; k<3; k++)
                        conv1[k] += colPixel[k];
                    nc += colPixel[0]+colPixel[1]+colPixel[2];


                    in.readRawData((char*)&valc[0], nimg);
                    vecc.assign(&valc[0], &valc[0]+nimg);
                    for(int k=0; k<nimg; k++)
                        conv2[k] += vecc[k];
                }
                for(int k=0; k<=size[0]-wid; k++)
                    in.readRawData((char*)&valc[0], nimg);
                vecc.assign(&valc[0], &valc[0]+nimg);
            }



        for(int k=0; k<nimg; k++)
            //qDebug() << conv2[k];
            //qDebug() << conv1[0] << " " << conv1[1] << " " << conv1[2];
            //qDebug() << nc;

            for(int k=0; k<3; k++)
                col[k] = conv1[k]/(wid*hei);
        nc = nc/(wid*hei);
        for(int k=0; k<nimg; k++)
            arrayInterpolation(k, 2) = conv2[k]/(wid*hei);

    }

    rbfsetpoints(modelInterpolation, arrayInterpolation);
    float rr = ui->rbfSpinBox->value();
    rbfsetalgomultilayer(modelInterpolation, rr, 1, 10e-3);

    rbfbuildmodel(modelInterpolation, repInterpolation);

    float avgi=0, zeroval=0, ni=0, mini=66666, maxi=0;

    for(int i=0; i<=100; i++)
        for(int j=0; j<=100; j++)
        {

            float lx = 0.02*(i-50);
            float ly = 0.02*(j-50);
            //float lz = (1-lx*lx-ly*ly);

            //qDebug() << lx << " " << ly  << "\n";

            if (lx*lx+ly*ly <= 1) {

                double res = rbfcalc2(modelInterpolation, lx, ly);
                if (res > maxi) maxi=res;
                if(res<mini) mini=res;
                avgi=avgi+res;
                if(i==50 && j==50)  zeroval=res;
                ni=ni+1;
                if(type==2)
                    for (int k=0; k<3; k++)
                        val[k] = 3*col[k]*res/(nc*256);
                if(type==1)
                    for (int k=0; k<3; k++)
                        val[k] = 3*col[k]*res/(nc);
            }

            else {
                val[0]=255;
                val[1]=255;
                val[2]=255;
            }

            for (int i2 = i*scale; i2 < i*scale+scale; i2++)
                for (int j2 = j*scale; j2 < j*scale+scale; j2++)
                    shownImg.at<Vec3b>(j2,i2) = val;
        }

    avgi=avgi/ni;
    ui->msgBox->setText("Plotted AP. Values at (0,0): " + QString::number(zeroval) + " average: " + QString::number(avgi));
    ui->msgBox->append("Min: " + QString::number(mini) + " max: " + QString::number(maxi));
    cv::flip(shownImg,shownImg,0);
    imshow (windowName, shownImg);

    filed.close();
}

void apTool::on_lxSpinBox_valueChanged(double i)
{

    int ly= (ui->lySpinBox->value())*(ui->lWidget->height()/(double)2);
    int lx= (ui->lxSpinBox->value())*(ui->lWidget->width()/(double)2);

    if (ly<=0) {
        ly = -ly;
        ly += ui->lWidget->height()/(double)2;
    }
    else {
        ly = ui->lWidget->height()/(double)2-ly;
    }

    if (lx>=0) {
        lx += ui->lWidget->width()/(double)2;
    }
    else {
        lx = ui->lWidget->width()/(double)2+lx;
    }

    QPixmap image = drawLWidget();
    QPainter painter(&image);

    painter.setBrush(QBrush(Qt::yellow));
    QRectF rectangle(lx-(5), ly-(5), 10.0, 10.0);
    painter.drawEllipse(rectangle);
    ui->lWidget->setPixmap(image);
}

void apTool::on_lySpinBox_valueChanged(double i)
{

    int ly= (ui->lySpinBox->value())*(ui->lWidget->height()/(double)2);
    int lx= (ui->lxSpinBox->value())*(ui->lWidget->width()/(double)2);

    if (ly<=0) {
        ly = -ly;
        ly += ui->lWidget->height()/(double)2;
    }
    else {
        ly = ui->lWidget->height()/(double)2-ly;
    }

    if (lx>=0) {
        lx += ui->lWidget->width()/(double)2;
    }
    else {
        lx = ui->lWidget->width()/(double)2+lx;
    }

    QPixmap image = drawLWidget();
    QPainter painter(&image);
    painter.setBrush(QBrush(Qt::yellow));
    QRectF rectangle(lx-(5), ly-(5), 10.0, 10.0);
    painter.drawEllipse(rectangle);
    ui->lWidget->setPixmap(image);
}

QPixmap apTool::drawLWidget() {

    QPixmap image(151, 151);
    QRectF rectangle;
    rectangle.setHeight(150.0);
    rectangle.setWidth(150.0);
    image.fill(Qt::lightGray);
    QPainter painter(&image);
    QColor color1(200, 120, 0);
    QColor color2(100, 50, 0);
    QRadialGradient gradient(rectangle.center(), 75);
    gradient.setColorAt(0, color1);
    gradient.setColorAt(1, color2);
    QPen pen(color2);
    painter.setPen(pen);
    QBrush brush(gradient);
    painter.setBrush(brush);
    painter.drawEllipse(rectangle);
    return image;

}

void apTool::on_fullTestBut_clicked()
{
    for(int i=0; i<54; i++){


        ui->robustMenu->setCurrentIndex(3);
        ui->viewBox->setCurrentIndex(4);
        ui->spinBox->setValue(i);
        ui->fitterMenu->setCurrentIndex(0);
        ui->processButton->click();
        ui->fitterMenu->setCurrentIndex(3);
        ui->processButton->click();

        ui->rbfSpinBox->setValue(0.3);
        ui->showButton->click();
        ui->rbfSpinBox->setValue(0.6);
        ui->showButton->click();

    }
}

void apTool::on_pushButton_8_clicked()
{
    QString fileName;

    fileName = QFileDialog::getOpenFileName(this,
                                            tr("Open grid file"));
    ui->gridEdit->setText(fileName);
}

void apTool::on_resampleButton_clicked()
{
    QString filename = ui->gridEdit->text();
    QFile file(filename);
    QTextStream textStream(&file);
    QString line="aa";
    QStringList parts;

    qDebug() << filename;

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    ui->viewBox->setCurrentIndex(2);
    int i=1;
    while(!(line.isNull())){
        line = textStream.readLine();
        parts = line.split(" ");
        double dire[3];
        double nofa=sqrt(parts[0].toFloat()*parts[0].toFloat()+parts[1].toFloat()*parts[1].toFloat()+parts[2].toFloat()*parts[2].toFloat());
        dire[0]=parts[0].toFloat()/nofa;
        dire[1]=parts[1].toFloat()/nofa;
        dire[2]=parts[2].toFloat()/nofa;

        this->ui->lxSpinBox->setValue(dire[0]);
        this->ui->lySpinBox->setValue(dire[1]);
         ui->showButton->click();


         QFile::rename("result.png","grid_" + QString::number(i) + ".png" );

         i++;
    }


}

void apTool::on_viewBox_activated(const QString &arg1)
{

}
