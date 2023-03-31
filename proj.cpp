#include "exo-vtk-include.h"
#include "config.h"
#include "mpi.h"

#define BIG
//#define SMALL

//*
//#define FICHIER MY_MESHES_PATH "/Frog_CHAR_X_256_Y_256_Z_44.raw"
//#define CHAR
//int startexploreval = 60;
//int endexploreval = 80;//*/


//#define FICHIER MY_MESHES_PATH "/Mystere5_SHORT_X_2048_Y_2048_Z_756.raw"
//
// int gridSize = 2048;
// int YgridSize = 2048;
// int ZgridSize = 756;
//
// #define SHORT
//
//int startexploreval=100;
//int endexploreval=65000;//


//voiture
//#define FICHIER MY_MESHES_PATH "/Mystere6_CHAR_X_1118_Y_2046_Z_694.raw"
//
//int gridSize = 1118;
//int YgridSize = 2046;
//int ZgridSize = 694;
//
//#define CHAR
//
//int startexploreval=13;
//int endexploreval=16;//


//tirelire
//#define FICHIER MY_MESHES_PATH "/Mystere1_SHORT_X_512_Y_512_Z_134.raw"
//int gridSize = 512;
//int YgridSize = 512;
//int ZgridSize = 134;
//
// #define SHORT
//
// int startexploreval=20000;
// int endexploreval=25000;//


//arbres
//#define FICHIER MY_MESHES_PATH "/Mystere2_SHORT_X_512_Y_400_Z_512.raw"
// #define SHORT
// int startexploreval=25000;
// int endexploreval=65000;//

//homard
//#define FICHIER MY_MESHES_PATH "/Mystere10_CHAR_X_1204_Y_1296_Z_224.raw"
//#define CHAR
//int startexploreval=50;
//int endexploreval=255;//

//poisson
//#define FICHIER MY_MESHES_PATH "/Mystere11_SHORT_X_512_Y_512_Z_1024.raw"
//#define SHORT
//int startexploreval=47000;
//int endexploreval=65000;//

//dent
//#define FICHIER MY_MESHES_PATH "/Mystere4_SHORT_X_512_Y_512_Z_322.raw"
//int gridSize = 512;
//int YgridSize = 512;
//int ZgridSize = 322;
//
// #define SHORT
//
// int startexploreval=50000;
// int endexploreval=65000; //

//tete
#define FICHIER MY_MESHES_PATH "/Mystere8_CHAR_X_2048_Y_2048_Z_2048.raw"
#define CHAR
int startexploreval=100;
int endexploreval=255; //


const char *location = FICHIER;

int gridSize;
int YgridSize;
int ZgridSize;

int winSize = 500;

int numPasses = 10;
int nbimages = 10;


const char *prefix = "";


int passNum = 0;

using std::cerr;
using std::endl;

// Function prototypes
vtkRectilinearGrid  *ReadGrid(int zStart, int zEnd);

void WriteImage(const char *name, const float *rgba, int width, int height);

bool ComposeImageZbuffer(float *rgba_out, float *zbuffer, int image_width, int image_height);

int strToInt(std::string str) {
    int result = 0;
    for (int i = 0; i < str.length(); i++) {
        result = result * 10 + (str[i] - '0');
    }
    return result;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    srand(time(NULL));
    std::string filename = location;
    gridSize = strToInt(
            filename.substr(filename.find("X_") + 2, filename.find("_Y_") - filename.find("X_") - 2));
    YgridSize = strToInt(
            filename.substr(filename.find("Y_") + 2, filename.find("_Z_") - filename.find("Y_") - 2));
    ZgridSize = strToInt(
            filename.substr(filename.find("Z_") + 2, filename.find(".raw") - filename.find("Z_") - 2));
    int npixels = winSize * winSize;
    vtkRectilinearGrid *reader = NULL;
    vtkLookupTable *lut = vtkLookupTable::New();
    lut->SetHueRange(0.1, 0.0);
    lut->SetSaturationRange(0.0, 1.0);
    lut->SetValueRange(1.0, 255.0);
    lut->SetNumberOfColors(100);
    lut->Build();
    vtkRenderer *ren = vtkRenderer::New();
    double bounds[6] = {0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001, 0.00001, 1 - 0.00001};
    ren->ResetCamera(bounds);
    vtkRenderWindow *renwin = vtkRenderWindow::New();
    renwin->SetSize(winSize, winSize);
    renwin->AddRenderer(ren);
    int zSize = ZgridSize/size;
    int zStart = rank * zSize;
    int zEnd = zStart + zSize;
    if (rank == size - 1) {
        zEnd = ZgridSize;
    }
    numPasses = 4;
    int passSize = (zEnd - zStart) / numPasses;
    int localStart = zStart;
    int localEnd = localStart + passSize;
    reader = ReadGrid(localStart, localEnd);

    vtkContourFilter *cf = vtkContourFilter::New();
    cf->SetNumberOfContours(1);
    int valcont = startexploreval;
    cf->SetValue(0, valcont);
    cf->SetInputData(reader);
    int maxsize = std::max(gridSize, std::max(YgridSize, ZgridSize));
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Scale(gridSize / (float) maxsize, YgridSize / (float) maxsize, ZgridSize / (float) maxsize);
    vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
    transformFilter->SetInputConnection(cf->GetOutputPort());
    transformFilter->SetTransform(transform);
    vtkDataSetMapper *mapper = vtkDataSetMapper::New();
    mapper->SetInputConnection(transformFilter->GetOutputPort());
    vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);
    mapper->SetScalarRange(startexploreval, endexploreval);
    mapper->SetLookupTable(lut);
    ren->AddActor(actor);
    ren->SetViewport(0, 0, 1, 1);
    renwin->SetOffScreenRendering(1);
    vtkCamera *cam = ren->GetActiveCamera();
    renwin->Render();


//    cam->Elevation(90); // si vous avez une valeure connu "à l'avance", mettez la avant la for
    // changer ici le nombre de rendu différent
    for(int az = 0; az<1; az++){
//        cf->SetValue(0, 30*az); // on la coupe ...
//        mapper->Update(); // sans oublier d'update le mapper en cas de modification de la coupe
        cam->Azimuth(36); // vous pouvez modifier les valeurs d'angles ici
//        cam->Elevation(36); // pareil pour l'élévation
        renwin->Render();
        float* rgba = new float[npixels * 4];
        float* zbuffer = new float[npixels];
        for(int px =0; px<npixels; px++){
            rgba[4*px] = 0;
            rgba[4*px+1] = 0;
            rgba[4*px+2] = 0;
            rgba[4*px+3] = 0;
            zbuffer[px] = 1.0;
        }
        for(int pass=0; pass<numPasses; pass++) {
            localStart = zStart + pass * passSize;
            localEnd = localStart + passSize;
            if(pass == numPasses - 1){
                localEnd = zEnd;
            }
            reader->Delete();
            reader = ReadGrid(localStart, localEnd);
            cf->SetInputData(reader);
            mapper->Update();
            renwin->Render();
            float* rgba2 = renwin->GetRGBAPixelData(0, 0, winSize - 1, winSize - 1, 0);
            float* zbuffer2 = renwin->GetZbufferData(0, 0, winSize - 1, winSize - 1);
            for(int i=0; i<npixels; i++){
                if(zbuffer[i] > zbuffer2[i]){
                    rgba[4*i] = rgba2[4*i];
                    rgba[4*i+1] = rgba2[4*i+1];
                    rgba[4*i+2] = rgba2[4*i+2];
                    rgba[4*i+3] = rgba2[4*i+3];
                    zbuffer[i] = zbuffer2[i];
                }
            }
            delete [] rgba2;
            delete [] zbuffer2;
        }
//        std::string name = std::to_string(rank)+"_az"+std::to_string(az)+".png";
//        WriteImage(name.c_str(), rgba, winSize, winSize);

        float* allrgba = nullptr;
        float* allzbuffer = nullptr;
        if(rank==0){
            allrgba = new float[4 * npixels * size];
            allzbuffer = new float[npixels *size ];
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(rgba, 4 * npixels, MPI_FLOAT, allrgba, 4 * npixels, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Gather(zbuffer, npixels, MPI_FLOAT, allzbuffer, npixels, MPI_FLOAT, 0, MPI_COMM_WORLD);
        delete [] rgba;
        delete [] zbuffer;

        if(rank==0){
            auto *finalrgba = new float[npixels * 4];
            auto *finalzbuffer = new float[npixels];
            for (int i = 0; i < npixels; i++) {
                finalzbuffer[i] = 1.0;
                finalrgba[i * 4] = 0.0;
                finalrgba[i * 4 + 1] = 0.0;
                finalrgba[i * 4 + 2] = 0.0;
                finalrgba[i * 4 + 3] = 1.0;
            }
            for (int i = 0; i < size; i++) {
                for (int index = 0; index < npixels; index++) {
                    int index2 = i * npixels + index;
                    if (allzbuffer[index2] < finalzbuffer[index]) {
                        finalzbuffer[index] = allzbuffer[index2];
                        finalrgba[index * 4] = allrgba[index2 * 4];
                        finalrgba[index * 4 + 1] = allrgba[index2 * 4 + 1];
                        finalrgba[index * 4 + 2] = allrgba[index2 * 4 + 2];
                        finalrgba[index * 4 + 3] = allrgba[index2 * 4 + 3];
                    }

                }
            }
            std::string finalfile = "final"+std::to_string(az)+".png";
            WriteImage(finalfile.c_str(), finalrgba, winSize, winSize);
        }
    }

    reader->Delete();
    mapper->Delete();
    cf->Delete();
    ren->RemoveActor(actor);
    actor->Delete();
    ren->Delete();
    renwin->Delete();
    MPI_Finalize();
}



// You should not need to modify these routines.
vtkRectilinearGrid *
ReadGrid(int zStart, int zEnd) {
    int i;

    /*   if (zStart < 0 || zEnd < 0 || zStart >= gridSize || zEnd >= gridSize || zStart > zEnd)
     {
     cerr << prefix << "Invalid range: " << zStart << "-" << zEnd << endl;
     return NULL;
     }
     */
    std::ifstream ifile(location);
    if (ifile.fail()) {
        cerr << prefix << "Unable to open file: " << location << "!!" << endl;
        throw std::runtime_error("can't find the file!! Check the name and the path of this file? ");
    }

//    cerr << prefix << "Reading from " << zStart << " to " << zEnd << endl;

    vtkRectilinearGrid *rg = vtkRectilinearGrid::New();
    vtkFloatArray *X = vtkFloatArray::New();
    X->SetNumberOfTuples(gridSize);
    for (i = 0; i < gridSize; i++)
        X->SetTuple1(i, i / (gridSize - 1.0));
    rg->SetXCoordinates(X);
    X->Delete();
    vtkFloatArray *Y = vtkFloatArray::New();
    Y->SetNumberOfTuples(YgridSize);
    for (i = 0; i < YgridSize; i++)
        Y->SetTuple1(i, i / (YgridSize - 1.0));
    rg->SetYCoordinates(Y);
    Y->Delete();
    vtkFloatArray *Z = vtkFloatArray::New();
    int numSlicesToRead = zEnd - zStart + 1;
    Z->SetNumberOfTuples(numSlicesToRead);
    for (i = zStart; i <= zEnd; i++)
        Z->SetTuple1(i - zStart, i / (ZgridSize - 1.0));
    rg->SetZCoordinates(Z);
    Z->Delete();

    rg->SetDimensions(gridSize, YgridSize, numSlicesToRead);

    unsigned int valuesPerSlice = gridSize * YgridSize;

#if defined(SHORT)
    unsigned int bytesPerSlice   = sizeof(short)*valuesPerSlice;

#elif defined(CHAR)
    unsigned int bytesPerSlice = sizeof(char) * valuesPerSlice;

#elif  defined(FLOAT)
    unsigned int bytesPerSlice   = sizeof(float)*valuesPerSlice;

#else
#error Unsupported choice setting
#endif


#if defined(SMALL)
    unsigned int offset = (unsigned int) zStart * (unsigned int) bytesPerSlice;
    unsigned int bytesToRead = bytesPerSlice * numSlicesToRead;
    unsigned int valuesToRead = valuesPerSlice * numSlicesToRead;
#elif defined(BIG)
    unsigned long long offset          = (unsigned long long)zStart * bytesPerSlice;
    unsigned long long  bytesToRead     = (unsigned long long )bytesPerSlice*numSlicesToRead;
    unsigned int valuesToRead    = (unsigned int )valuesPerSlice*numSlicesToRead;
#else
#error Unsupported choice setting
#endif


#if defined(SHORT)
    vtkUnsignedShortArray *scalars = vtkUnsignedShortArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    unsigned short *arr = scalars->GetPointer(0);

#elif defined(CHAR)
    vtkUnsignedCharArray *scalars = vtkUnsignedCharArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    unsigned char *arr = scalars->GetPointer(0);

#elif  defined(FLOAT)
    vtkFloatArray *scalars = vtkFloatArray::New();
    scalars->SetNumberOfTuples(valuesToRead);
    float *arr = scalars->GetPointer(0);
#else
#error Unsupported choice setting
#endif


    ifile.seekg(offset, ios::beg);
    ifile.read((char *) arr, bytesToRead);
    ifile.close();

    scalars->SetName("entropy");
    rg->GetPointData()->AddArray(scalars);
    scalars->Delete();

    vtkFloatArray *pr = vtkFloatArray::New();
    pr->SetNumberOfTuples(valuesToRead);
#if defined(SMALL)
    for (unsigned int i = 0; i < valuesToRead; i++)
#elif defined(BIG)
        for (unsigned long long  i = 0 ; i < valuesToRead ; i++)
#endif
            pr->SetTuple1(i, passNum);

    pr->SetName("pass_num");
    rg->GetPointData()->AddArray(pr);
    pr->Delete();

    rg->GetPointData()->SetActiveScalars("entropy");

//    cerr << prefix << " Done reading" << endl;
    return rg;
}
void
WriteImage(const char *name, const float *rgba, int width, int height) {
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
#if VTK_MAJOR_VERSION <= 5
    img->SetNumberOfScalarComponents(3);
    img->SetScalarTypeToUnsignedChar();
#else
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
#endif

    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++) {
            unsigned char *ptr = (unsigned char *) img->GetScalarPointer(i, j, 0);
            int idx = j * width + i;
            ptr[0] = (unsigned char) (255 * rgba[4 * idx + 0]);
            ptr[1] = (unsigned char) (255 * rgba[4 * idx + 1]);
            ptr[2] = (unsigned char) (255 * rgba[4 * idx + 2]);
        }


    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(name);
    writer->Write();

    img->Delete();
    writer->Delete();
}