#include<cmath>
#include<vector>
#include<iostream>
#include<string>


struct Point {
    std::vector<double> coords;
    double cost;
    Point *parent;

    Point(std::vector<double> coords, Point *parent=NULL, double cost=0): coords(coords), parent(parent), cost(cost) {
    }

    Point(int n, double cost=0): coords(n,0), cost(cost){
    }

    Point() {
    }

    Point(const Point &point): coords(point.coords), parent(point.parent), cost(point.cost) {
    
    }

    double sqDistance(std::vector<double> &other) {
        double dist = 0;
        for(int i=0; i<coords.size(); i++){
            dist += (coords[i]-other[i])*(coords[i]-other[i]);
        }
        return dist;
    }

    double sqDistance(Point& other) {
        return sqDistance(other.coords);
    }

};

class KDTree {

    static constexpr double eps = 1e-5;
    unsigned long long treeSize;
    int dim;

    struct Node {
        Point point;
        Node *left;
        Node *right;
        Node(Point x) : point(x), left(NULL), right(NULL) {}
        Point* getPointAddress(){
            return &point;
        }
    } *root;

    Node* newNode(Point &point) {
        Node *node = new Node(point);
        treeSize++;
        return node;
    }

    void insertRec(Node *root, Point &point, unsigned depth, Point* &final){
    
        // Calculate current dimension (cd) of comparison
        unsigned cd = depth % dim;
    
        // Compare the new point with root on current dimension 'cd'
        // and decide the left or right subtree
        if (point.coords[cd] < (root->point.coords[cd])){
            if(root->left == NULL){
                root->left  = newNode(point);
                final = root->left->getPointAddress();
            }
            else
                insertRec(root->left, point, depth + 1, final);
        }
        else{
            if(root->right == NULL){
                root->right = newNode(point);
                final = root->right->getPointAddress();
            }
            else
                insertRec(root->right, point, depth + 1, final);
        }
    }

    Point* searchRec(Node* root, Point &point, unsigned depth){
        // Base cases
        if (root == NULL)
            return NULL;
        if (root->point.sqDistance(point) < eps)
            return root->getPointAddress();
    
        // Current dimension is computed using current depth and total
        // dimensions (k)
        unsigned cd = depth % dim;
    
        // Compare point with root with respect to cd (Current dimension)
        if (point.coords[cd] < root->point.coords[cd])
            return searchRec(root->left, point, depth + 1);
    
        return searchRec(root->right, point, depth + 1);
    }

    void nearestNeighborRec(std::vector<double> &point, Node *root, Point* &guess, double &minDistSq, unsigned depth=0){
        if(root == NULL)
            return;
        
        unsigned cd = depth % dim;

        double dist = root->point.sqDistance(point);
        if(dist < minDistSq){
            minDistSq = dist;
            guess = root->getPointAddress();
        }

        double dx = fabs(point[cd] - root->point.coords[cd]);
        dx = dx*dx;

        if(point[cd] < root->point.coords[cd]){
            nearestNeighborRec(point, root->left, guess, minDistSq, depth+1);
            if(dx <= minDistSq)
                nearestNeighborRec(point, root->right, guess, minDistSq, depth+1);
        } else {
            nearestNeighborRec(point, root->right, guess, minDistSq, depth+1);
            if(dx <= minDistSq)
                nearestNeighborRec(point, root->left, guess, minDistSq, depth+1);
        }
    }

    void searchRadiusRec(std::vector<double> &point, Node *root, double radius, std::vector<Point*> &result, unsigned depth=0){
        if(root == NULL)
            return;
        
        unsigned cd = depth % dim;
        double dist = root->point.sqDistance(point);
        if(dist <= radius*radius){
            result.push_back(root->getPointAddress());
        }

        double dx = fabs(point[cd] - root->point.coords[cd]);
        if(point[cd] < root->point.coords[cd]){
            searchRadiusRec(point, root->left, radius, result, depth+1);
            if(dx <= radius)
                searchRadiusRec(point, root->right, radius, result, depth+1);
        } else {
            searchRadiusRec(point, root->right, radius, result, depth+1);
            if(dx <= radius)
                searchRadiusRec(point, root->left, radius, result, depth+1);
        }
    }

    void exportTreeRec(std::ostream &out, Node *root){
        if(root == NULL)
            return;
        for(auto &x: root->point.parent->coords)
            out<<x<<" ";
        for(auto &x: root->point.coords)
            out<<x<<" ";
        out<<std::endl;
        exportTreeRec(out, root->left);
        exportTreeRec(out, root->right);
    }

    public:

    KDTree(int dim): dim(dim), root(NULL), treeSize(0) {
    }

    Point* insert(Point point) {
        if (root==NULL) {
            root = newNode(point);
            return &root->point;
        }

        Point* final = NULL;
        insertRec(root, point, 0, final);
        return final;
    }

    Point* search(Point point) {
        return searchRec(root, point, 0);
    }

    Point* nearestNeighbor(std::vector<double> point, double* minDist=NULL){
        Point *guess = NULL;
        double _minDist = 1e9;

        nearestNeighborRec(point, root, guess, _minDist);

        if(minDist)
            *minDist = _minDist;
        
        return guess;
    }

    std::vector<Point*> searchRadius(std::vector<double> point, double radius){
        std::vector<Point*> result;
        searchRadiusRec(point, root, radius, result);
        return result;
    }

    Node* getRoot(){
        return root;
    }

    unsigned long long size(){
        return treeSize;
    }

    void exportTree(std::string out){
        std::ofstream fout(out);
        exportTreeRec(fout, root->left);
        exportTreeRec(fout, root->right);
    }
};

/*

int main(){
    // Create a KDTree object with dimension 2
    KDTree tree(2);

    // Create some points
    std::vector<std::vector<double>> points = {{3, 6}, {17, 15}, {13, 15}, {6, 12}, {9, 1}, {2, 7}, {10, 19}};

    // Insert the points into the KDTree
    for(auto &point : points){
        tree.insert(Point(point));
    }

    // Search for a point
    std::vector<double> searchPoint = {3, 6};
    Point* result = tree.search(Point(searchPoint));
    if(result != NULL){
        std::cout << "Found point at (" << result->coords[0] << ", " << result->coords[1] << ")" << std::endl;
    } else {
        std::cout << "Point not found" << std::endl;
    }

    // Find the nearest neighbor
    std::vector<double> queryPoint = {2.0, 4.0};
    double minDist;
    Point* nearestNeighbor = tree.nearestNeighbor(queryPoint, &minDist);
    if(nearestNeighbor != NULL){
        std::cout << "Nearest neighbor: (" << nearestNeighbor->coords[0] << ", " << nearestNeighbor->coords[1] << ")" << std::endl;
        std::cout << "Distance: " << minDist << std::endl;
    } else {
        std::cout << "No nearest neighbor found" << std::endl;
    }

    // Search within a radius
    std::vector<double> centerPoint = {12.0, 16.0};
    double radius = 6.0;
    std::vector<Point*> pointsWithinRadius = tree.searchRadius(centerPoint, radius);
    std::cout << "Points within radius: " << std::endl;
    for(auto point : pointsWithinRadius){
        std::cout << "(" << point->coords[0] << ", " << point->coords[1] << ")" << std::endl;
    }

    return 0;
}

*/