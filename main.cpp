#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageSelector.h>
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include "DGtal/io/Color.h"
#include <DGtal/topology/DigitalTopology.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/io/readers/PGMReader.h>
#include <DGtal/io/writers/PGMWriter.h>
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/geometry/curves/estimation/DSSLengthEstimator.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"

using namespace std;
using namespace DGtal;
using namespace Z2i;

typedef Object<DT4_8, DigitalSet> ObjectType; // Digital object type

template<class T>
std::vector<Z2i::SCell> getBoundary(T & object)
{
    //Khalimsky space
    typedef KhalimskySpaceND< 2, int > KSpace;
    KSpace K;
    // we need to add a margine to prevent situations such that an object touch the bourder of the domain
    K.init( object.domain().lowerBound() - Point(1,1),
                   object.domain().upperBound() + Point(1,1), true );
    
    // 1) Call Surfaces::findABel() to find a cell which belongs to the border
    //Extract a boundary cell
    Z2i::SCell aCell = Surfaces<Z2i::KSpace>::findABel(K, object.pointSet(),10000);
 
    std::vector<Z2i::SCell> vectBdrySCell;
    // 2) Call Surfece::track2DBoundaryPoints to extract the boundary of the object
    SurfelAdjacency<2> SAdj( true );
    Surfaces<Z2i::KSpace>::track2DBoundary( vectBdrySCell, K, SAdj, object.pointSet(), aCell );

    return vectBdrySCell;
}


template<class T>
void sendToBoard(Board2D & board, T & p_Object, DGtal::Color p_Color,int numberOfImage) {// show Digital Object
    board << CustomStyle(p_Object.className(), new DGtal::CustomFillColor(p_Color));
    board << p_Object;
    std::ostringstream path;
    path << "Digital object of rice grain from image number (" << numberOfImage << ").eps";
    board.saveEPS(path.str().c_str());
}

// Calculate the area of a polygon using the shoelace formula
double calculatePolygonArea(const std::vector<LibBoard::Point>& polygonVertices) {
    int n = polygonVertices.size(); 
    double area = 0.0; 

    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n; 

        // Apply the shoelace formula to calculate the area
        area += (polygonVertices[i].x * polygonVertices[j].y) - (polygonVertices[j].x * polygonVertices[i].y);
    }

    area = std::abs(area) / 2.0;
    return area; 
}

// Calculate the perimeter of a polygon by summing the distances between its vertices
double calculatePolygonPerimeter(const std::vector<LibBoard::Point>& polygonVertices) {
    double perimeter = 0.0; 
    int n = polygonVertices.size(); 

    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;

        // Calculate the distance between two adjacent vertices and add it to the perimeter
        perimeter += (polygonVertices[i] - polygonVertices[j]).norm();
    }

    return perimeter; 
}

// DSS takes a Curve object and returns a vector of points
std::vector<LibBoard::Point> DSS(Curve C, int numberOfImage) {

    typedef Z2i::Curve::PointsRange Range;
    Range R = C.getPointsRange(); 

    std::vector<Z2i::Point> directions = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};

    // Initialize the Freeman chain code string
    std::string freemanChainCode;

    Z2i::Curve::PointsRange::ConstIterator it = R.begin();

    Z2i::Point begin = *it; 
    Z2i::Point prevPoint = *it; 
    ++it;

    // Iterate through the points in the range
    for (; it != R.end(); ++it) {
        Z2i::Point currentPoint = *it;
        Z2i::Point displacement = currentPoint - prevPoint;

        int direction = -1;
        for (int i = 0; i < directions.size(); ++i) {
            if (directions[i] == displacement) {
                direction = i;
                break;
            }
        }

        if (direction != -1) {
            freemanChainCode += std::to_string(direction);
        }

        prevPoint = currentPoint;
    }

    // Calculate the direction for the last displacement to close the circle
    Z2i::Point d = begin - prevPoint;
    int direction = -1;
    for (int i = 0; i < directions.size(); ++i) {
        if (directions[i] == d) {
            direction = i;
            break;
        }
    }

    if (direction != -1) {
        freemanChainCode += std::to_string(direction);
    }


    
    typedef FreemanChain<int> Contour4;
    typedef ArithmeticalDSSComputer<Contour4::ConstIterator, int, 4> DSS4;
    typedef GreedySegmentation<DSS4> Decomposition4;

    std::stringstream ss(stringstream::in | stringstream::out);
    ss << "0 0 " + freemanChainCode << std::endl;

    typedef FreemanChain<int> Contour4;
    Contour4 theContour(ss);

    Decomposition4 theDecomposition(theContour.begin(), theContour.end(), DSS4());

    
    Board2D dBoard;
    dBoard << SetMode("PointVector", "Grid");

    // Create a vector to store LibBoard::Point points
    std::vector<LibBoard::Point> points;
    LibBoard::Point firstPoint(100000000, 100000000);

    
    for (Decomposition4::SegmentComputerIterator it = theDecomposition.begin(), itEnd = theDecomposition.end(); it != itEnd; ++it) {
        if (firstPoint.x == 100000000) {
            firstPoint.x = (*it).Ul()[0];
            firstPoint.y = (*it).Ul()[1];
        }

        points.push_back(LibBoard::Point((*it).Ul()[0], (*it).Ul()[1]));

        // dBoard << SetMode("ArithmeticalDSS", "Points") << it->primitive();
        // dBoard << SetMode("ArithmeticalDSS", "BoundingBox") << CustomStyle("ArithmeticalDSS/BoundingBox", new CustomPenColor(Color::Blue)) << it->primitive();
    }

    // Add the first point again to close the polyline
    points.push_back(LibBoard::Point(firstPoint.x, firstPoint.y));


    // std::ostringstream path;
    // path << "DSS of Rice grain from image number (" << numberOfImage << ").eps";
    // dBoard.saveEPS(path.str().c_str());


    return points;
}

std::vector< ObjectType > Eliminate(std::vector< ObjectType > objs,DGtal::Z2i::Domain domain) {
        // Create a vector to store the indices of objects to be deleted.
    std::vector<int> positionObjectsToDelete;

    for (int i = 0; i < objs.size(); ++i) {
        ObjectType& obj = objs[i]; 

        for (const Z2i::Point& point : obj.pointSet()) {
            // Check if the point is close to the border.
            if (point[0] <= 0|| point[0] >= domain.upperBound()[0]||
                point[1] <= 0|| point[1] >= domain.upperBound()[1]) {
              
                positionObjectsToDelete.push_back(i);
                
                break; 
            }
        }
    }

    // Remove marked objects from 'objs'.
    for (int i = positionObjectsToDelete.size() - 1; i >= 0; --i) {
        objs.erase(objs.begin() + positionObjectsToDelete[i]);
    }
    return objs;
}

int main(int argc, char** argv)
{
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings
    typedef ImageSelector<Domain, unsigned char >::Type Image; // type of image
    typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet; // Digital set type
    typedef Domain::ConstIterator DomainConstIterator; 


    /////////////////////////////////////// step 2 //////////////////////////////////////////////
    std::vector<std::string> imageFile = {
        "./RiceGrains/Rice_japonais_seg_bin.pgm",
        "./RiceGrains/Rice_camargue_seg_bin.pgm",
        "./RiceGrains/Rice_basmati_seg_bin.pgm"
        };
    int numberOfImage=0;
    for (const std::string &imageP : imageFile) {
        // read an image
    Image image = PGMReader<Image>::importPGM (imageP); 

    // 1) make a "digital set" of proper size
    Z2i::DigitalSet digitalSet (image.domain());
    // 2) populate a digital set from the image using SetFromImage::append()
    SetFromImage<Z2i::DigitalSet>::append<Image>(digitalSet, image, 1, 255);

    Board2D aBoard;
    aBoard << digitalSet;
    // aBoard << image.domain();
    std::ostringstream path;
    path << "Digital Set of image number (" << numberOfImage << ").eps";
    aBoard.saveEPS(path.str().c_str()); // save image 

    
    typedef SpaceND< 2,int > Z2;
    typedef MetricAdjacency< Z2, 1 > Adj4;
    typedef MetricAdjacency< Z2, 2 > Adj8;
    typedef DigitalTopology< Adj4, Adj8 > DT4_8;
    
    Adj4 adj4;
    Adj8 adj8;
    DT4_8 dt4_8( adj4, adj8);

    // 3) Create a digital object from the digital set
    std::vector< ObjectType > objs;          // All connected components are going to be stored in it
    std::back_insert_iterator< std::vector< ObjectType > > inserter( objs ); // Iterator used to populated "objects".

    ObjectType objects( dt4_8, digitalSet);
    // 4) Set the adjacency pair and obtain the connected components using "writeComponents"
    objects.writeComponents( inserter );
    std::cout << "number of digital object of image number "<< numberOfImage << " : " << objs.size() << endl; // Right now size of "objects" is the number of conected components

    ////////////////// eliminate the grains whose entire bodies do not appear in the image frame///////////
    objs = Eliminate(objs,image.domain());
    std::cout << "number of digital object of image number after the elimination"<< numberOfImage << " : " << objs.size() << endl;
    ///////////'objs' now contains the objects that are not partially on the border.///////
  

    Board2D bBoard;
    sendToBoard(bBoard, objs[10], Color::Red, numberOfImage);   // send the connected component "objs[10]" to "bBoard"

/////////////////////////////////// Step3 /////////////////////////////////////////////////
    Board2D cBoard;
    std::vector<vector<Z2i::SCell>> vectBdry; //// store the boundry of all grains in the image
    std::vector<ObjectType>::iterator objectIterator;
    for ( objectIterator=objs.begin() ; objectIterator!= objs.end(); objectIterator++ ){    
        vectBdry.push_back(getBoundary(*objectIterator));
    }

    Board2D board;

    std::vector<Z2i::SCell>::iterator it;
    for ( it=vectBdry.at(10).begin() ; it != vectBdry.at(10).end(); it++ ){
        board<< CustomStyle((*it).className() ,new DGtal::CustomFillColor( Color::Red))<< *it;
    }

 
    std::ostringstream path1;
    path1 << "Digital Object Boundery of rice grain from image number (" << numberOfImage << ").eps";
    board.saveEPS(path1.str().c_str()); //// save image

///////////////////////////////////////// step 4 (Polygonization)///////////////////////////////////////////       
    std::vector<Z2i::Curve > Cvector ; //// store the curve of all grains in the image
    std::vector<vector<Z2i::SCell>>::iterator bdryIterator;
     for ( bdryIterator=vectBdry.begin() ; bdryIterator!= vectBdry.end(); bdryIterator++ ){
            Z2i::Curve boundaryCurve;
            boundaryCurve.initFromSCellsVector(*bdryIterator);
            Cvector.push_back(boundaryCurve);
    }


    std::vector<std::vector<LibBoard::Point >> vectorPoints; //// store the points of polygones from all the grains in the image

    std::vector<Z2i::Curve>::iterator pointIterator;

    for ( pointIterator=Cvector.begin() ; pointIterator!= Cvector.end(); pointIterator++ ){

        std::vector<LibBoard::Point > points = DSS(*pointIterator,numberOfImage); 
        vectorPoints.push_back(points);
    }
    Board2D eBoard;
    eBoard.setPenColor(Color::Red);

    // Draw the polyline of the polygone 
    eBoard.drawPolyline(vectorPoints.at(10));

    
    std::ostringstream path2;
    path2 << "DSS of Rice grain from image number (" << numberOfImage << ").eps";
    eBoard.saveEPS(path2.str().c_str()); // save image


///////////////////////////////// step 5 (Area) ////////////////////////////

///////////////////////////////////2-cells ////////////////////////////////////
std::cout << "Area of Digital Object (2-cells) of image : "<< numberOfImage<< std::endl;
string sA2=" ";
double moyA2;
for (const ObjectType& obj : objs) {
    DigitalSet componentSet = obj.pointSet();
    std::size_t area = componentSet.size();
    sA2 = sA2 + " , " + std::to_string(area);
    moyA2 = moyA2+area;
    
}
std::cout << "Area = " << sA2<< std::endl;
moyA2 = moyA2/objs.size();
std::cout << "moyenne of Area of Digital Object (2-cells) of image "<< numberOfImage<< " = " << moyA2 << std::endl;

///////////////////////////   Shoelace /////////////////////////////////////
double moySh = 0;
string sSh=" ";
std::cout << "Area of Digital Object (shoelace formula) of image : "<< numberOfImage<< std::endl;
for (std::vector<LibBoard::Point>& point : vectorPoints) {
double polygonArea = calculatePolygonArea(point);
moySh = moySh+polygonArea;
sSh = sSh + " , " + std::to_string(polygonArea); 
}
std::cout << "Area (shoelace)  =  "  << sSh << std::endl;
moySh = moySh/vectorPoints.size();
std::cout << "moyenne of Area of Digital Object (shoelace formula) of image "<< numberOfImage<< " = " << moySh << std::endl;



//////////////////////////////////setp 6 //////////////////////////////////////////


//////////////////////// perimetre 1-cell /////////////////////////////////////////
double moyPer1 = 0;
string sPer1=" ";
std::cout << "Perimeter of Digital Object (1-Cells) of image : "<< numberOfImage<< std::endl;
for (vector<Z2i::SCell>& bdry : vectBdry) {
    sPer1 = sPer1 + " , " + std::to_string(bdry.size()); 

moyPer1 = moyPer1+bdry.size();
}
std::cout << "Perimeter (1-cells)  =  "  << sPer1 << std::endl;
moyPer1 = moyPer1/vectorPoints.size();
std::cout << "moyenne of Perimeter of Digital Object (1-Cells) of image "<< numberOfImage<< " = " << moyPer1 << std::endl;

// ////////////////////////////////// as the perimeter of the polygon ////////////////////////////////


double moyPerSeg = 0;
string sPerSeg=" ";
std::cout << "Perimeter of Digital Object (polygone segment) of image : "<< numberOfImage<< std::endl;
for (std::vector<LibBoard::Point >& segments : vectorPoints) {
    double polygonPerimeter = calculatePolygonPerimeter(segments);
    sPerSeg = sPerSeg + " , " + std::to_string(polygonPerimeter); 

moyPerSeg = moyPerSeg+polygonPerimeter;
}
std::cout << "Perimeter (polygone segment)  =  "  << sPerSeg << std::endl;
moyPerSeg = moyPerSeg/vectorPoints.size();
std::cout << "moyenne of Perimeter of Digital Object (polygone segment) of image "<< numberOfImage<< " = " << moyPerSeg << std::endl;


// ////////////////////////////////step 7/////////////////////////////////////////////////////////////


double moyCir = 0;
string sCir = " ";
std::cout << "Circularity of digital objects  of image :"<< numberOfImage<< std::endl;
for(int i=0;i<objs.size();i++){
    double polygonArea = calculatePolygonArea(vectorPoints.at(i));
    double polygonPerimeter = calculatePolygonPerimeter(vectorPoints.at(i));
    double circularity = 4*M_PI*polygonArea/(polygonPerimeter*polygonPerimeter); 
    sCir = sCir + " , " + std::to_string(circularity);
    moyCir = moyCir + circularity;  
    
}
std::cout << "Circularity of digital objects = " << sCir << std::endl;
moyCir = moyCir / objs.size();
std::cout << "moyenne of Circularity of Digital Object of image "<< numberOfImage<< " = " << moyCir << std::endl;

numberOfImage++;
    }
    
  return 0;
}
