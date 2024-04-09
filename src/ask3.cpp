#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cmath>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

using namespace tbb;

struct Body;
class bhTreeNode;
double lengthOfSmallestGridContainingBody(bhTreeNode* root,  Body& b);

int REPETITIONS;
int numBodies;
double universeSize;
std::vector<Body> bodies;
bhTreeNode *TreeRoot;

// struct that contains the info of each body
class Body {
public:
    double x, y, Vx, Vy, mass, netForceX = 0, netForceY = 0;
    std::string name;

    void calculateNewPossition(){
        double ax = this->netForceX / this->mass;
        double ay = this->netForceY / this->mass;

        this->Vx = this->Vx + 1 * ax;
        this->Vy = this->Vy + 1 * ay;

        this->x = this->x + 1 * this->Vx;
        this->y = this->y + 1 * this->Vy;
    }


};

class bhTreeNode{
private:
    static constexpr double G = 6.67e-11; // Gravitational constant

public:
    double TotalMass;
    double centerMassX;
    double centerMassY;
    double minX, minY, maxX, maxY; // the array this node covers
    bool isLeaf;
    std::vector<bhTreeNode*> children; // one child for each quadrant


    explicit bhTreeNode(double minX = 0, double minY = 0, double maxX = 0, double maxY = 0)
            : TotalMass(0), centerMassX(0), centerMassY(0), minX(minX), minY(minY), maxX(maxX), maxY(maxY), isLeaf(false){}


    double length() const{return this->maxX - this->maxY;}

    bool isBodyInGrid(Body &b) const{
        if(b.x >= this->minX && b.x <= this->maxX && b.y >= this->minY && b.y <= this->maxY)
            return true;
        return false;
    }

    void netForce(Body &b){
            double F;
            double minLen = lengthOfSmallestGridContainingBody(TreeRoot, b);

            if((!(this->isBodyInGrid(b))) && ((this->centerMassX - b.x) > minLen) && ((this->centerMassY - b.y) > minLen)){
                // calculate with the aprox
                double dx = this->centerMassX - b.x;
                double dy = this->centerMassY - b.y;
                double r_squared = dx * dx + dy * dy;
                if(r_squared != 0){
                    F = (G * this->TotalMass * b.mass) / r_squared;

                    b.netForceX = (F * (b.x - this->centerMassX)) / sqrt(r_squared);
                    b.netForceY = (F * (b.y - this->centerMassY)) / sqrt(r_squared);
                }

            }else{

                //clasulate for each children
                for (bhTreeNode* child : children) {
                    child->netForce(b);
                }


                //not worth the overhead
    //                    parallel_for(blocked_range<size_t>(0, children.size()), [&](const blocked_range<size_t>& r) {
    //                        for(size_t i = r.begin(); i != r.end(); ++i) {
    //                            children[i]->netForce(b);
    //                        }
    //                    });
            }
    }

};


void readInput(const std::string& path){
    std::ifstream inputFile(path);
    if (!inputFile) {
        std::cerr << "Error opening file!" << std::endl;
        exit(1);
    }

    inputFile >> numBodies;
    inputFile >> universeSize;

    bodies.resize(numBodies);
    for (int i = 0; i < numBodies; ++i) {
        inputFile >> bodies[i].x >> bodies[i].y >> bodies[i].Vx >> bodies[i].Vy >> bodies[i].mass;
        std::getline(inputFile >> std::ws, bodies[i].name);
        // some name of the input finish with the \r char so we remove it.
        if (!bodies[i].name.empty() && bodies[i].name.back() == '\r') {
            bodies[i].name.pop_back();
        }
    }
}


// Function to create a quadtree from a vector of bodies
bhTreeNode* createTree(double minX, double minY, double maxX, double maxY) {
    // Create the root node
    auto* root = new bhTreeNode(minX, minY, maxX, maxY);

    // Find the bodies in our bounds
    long int bodyCount = 0;
    double sumOfxTimesMass = 0.0;
    double sumOfyTimesMass = 0.0;
    for (auto &body : bodies) {
        if (root->isBodyInGrid(body)) {
            bodyCount++;
            root->TotalMass += body.mass;
            sumOfxTimesMass += body.x * body.mass;
            sumOfyTimesMass += body.y * body.mass;

        }
    }

    // Calculate center of mass
    if (root->TotalMass > 0) {
        root->centerMassX = sumOfxTimesMass / root->TotalMass;
        root->centerMassY = sumOfyTimesMass / root->TotalMass;
    }

    // Base case
    if (bodyCount <= 1) {
        root->isLeaf = true;
        return root;
    }


    // Calculate midpoints
    double midX = (minX + maxX) / 2;
    double midY = (minY + maxY) / 2;

    // Create children (quadrants)

    // 1st quadrant (nw)
    root->children.push_back(createTree(minX, midY, midX, maxY));

    // 2nd quadrant (ne)
    root->children.push_back(createTree(midX, midY, maxX, maxY));

    // 3rd quadrant (sw)
    root->children.push_back(createTree(minX, minY, midX, midY));

    // 4th quadrant (se)
    root->children.push_back(createTree(midX, minY, maxX, midY));

    return root;
}


double lengthOfSmallestGridContainingBody(bhTreeNode* root,  Body& b) {
    if( root->isLeaf && root->isBodyInGrid(b) )return root->length();
    if(!(root->isBodyInGrid(b))) return std::numeric_limits<double>::max();

    // Recursively find the smallest grid containing the body in each child node
    //double minGridSize = std::numeric_limits<double>::max();
    double  minGridSize = root->length();
    for (bhTreeNode* child : root->children) {
        if (child->isBodyInGrid(b)) {
            double childSize = lengthOfSmallestGridContainingBody(child, b);
            minGridSize = std::min(minGridSize, childSize);
        }
    }
    return minGridSize;
}

// Function to print the quadtree
void printTree(bhTreeNode *root, int depth = 0) {
    if (root == nullptr) {
        return;
    }

    // Print node information
    std::cout << std::string(depth, '-') << " Node" << std::endl;
    std::cout << std::string(depth, '-') << " Total Mass: " << root->TotalMass << std::endl;
    std::cout << std::string(depth, '-') << " Bounds: (" << root->minX << ", " << root->minY << ") - (" << root->maxX << ", " << root->maxY << ")" << std::endl;

    // Print children
    if (!root->children.empty()) {
        std::cout << std::string(depth, '-') << " Children:" << std::endl;
        for (const auto &child : root->children) {
            printTree(child, depth + 1);
        }
    }
}

void runSimulation(){
    std::cout << "----------------------- Gravitational simulation using Barnes-Hut algorithm -----------------------" << std::endl;
    for(int n = 0; n < REPETITIONS; n++){
        //create Tree
        TreeRoot = createTree(-universeSize, -universeSize, universeSize, universeSize);

        // calculate forces and new position(they dont have race condition as each pressure of a bodied only read from the tree
        // and its self and only changes the data to its self and not the tree that the other bodies reed)
        parallel_for(

                blocked_range<size_t>(0, bodies.size()),

                [&](const blocked_range<size_t>& r) -> void {
                           for(size_t i = r.begin(); i != r.end(); i++ ) {
                               // Calculate net force acting on the body
                               TreeRoot->netForce(bodies[i]);
                               // Update position of the body according to the forces
                               bodies[i].calculateNewPossition();
                           }
                       }
                       
         );




        printf("\r %d/%d seconds simulated", n+1, REPETITIONS);
        fflush(stdout);
    }

    printf("\n\n");

    delete TreeRoot;
}

void printResults() {
    std::ofstream outputFile("output.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    outputFile << numBodies << std::endl;
    outputFile << universeSize << std::endl;

    for (const auto& body : bodies) {
        outputFile << body.x << " " << body.y << " " << body.Vx << " " << body.Vy << " " << body.mass << " " << body.name << std::endl;
    }

    outputFile.close();
}


int main (int argc, char **argv){
    std::string input;

    if(argc == 1){
        std::cout << "please give the path of the input file\n>";
        std::cin >> input;

        std::cout << "please give the number of repetitions to perform\n>";
        std::cin >> REPETITIONS;
    }
    else if(argc == 3) {
        input = argv[1]; REPETITIONS = atoi(argv[2]);
    }
    else {
        std::cerr << "invalid call\nPass no arguments to give the input from stdin or only two"
                     " (the paht of the input) and the number of repetitions\n";
        exit(-1);
    }

    readInput(input);

    runSimulation();

    printResults();

    return 0;
}
