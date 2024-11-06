#include<vector>
#include<random>
#include<cmath>
#include<iostream>
#include<fstream>

#include<omp.h>
#include<mpi.h>

#include "kdtree.cpp"
#include "collision.cpp"

#define N 

int myid;

class RRTStar {

    static std::mt19937 gen;
    std::uniform_real_distribution<> dis;
    omp_lock_t treeLock, sendLock;

    protected:
    
    KDTree nodes;
    CollisionChecker cc;
    std::vector<double> goalCoord;
    double stepSize;
    double goalRadius;
    double rewireRadius;
    int dof;

    std::vector<double> generateRandomPoint() {
        std::vector<double> point(dof,0);
        for(auto &x: point){
            x = dis(gen);
        }
        return point;
    }

    std::vector<double> moveStep(std::vector<double>& start, std::vector<double>& end){
        std::vector<double> res(end);

        for(int i=0; i<res.size(); i++)
            res[i] -= start[i];
        
        double mag = 0;
        for(auto x: res){
            mag += x*x;
        }

        if(mag<stepSize*stepSize)
            return end;
        
        mag = sqrt(mag);
        for(int i=0; i<res.size(); i++)
            res[i] = stepSize*res[i]/mag + start[i];

        return res;
    }

    static double eulerDistance(std::vector<double>& start, std::vector<double>& end){
        std::vector<double> res(end);

        for(int i=0; i<res.size(); i++)
            res[i] -= start[i];
        
        double mag = 0;
        for(auto x: res){
            mag += x*x;
        }
        return mag;
    }

    bool checkGoalState(std::vector<double>& point, std::vector<double>& otherCoord){

        bool result=true;

        for(int i=0; i<point.size(); i++){
            if(fabs(otherCoord[i]-point[i])>goalRadius)
                result = false;
        }
        
        return result;
    }

    Point* extend(){
        std::vector<double> rand_point = generateRandomPoint();

        Point *nn = nodes.nearestNeighbor(rand_point);
        std::vector<double> newPoint = moveStep(nn->coords, rand_point);
        
        if(checkCollision(newPoint)==true)
            return NULL;

        double cost = nn->cost + eulerDistance(nn->coords, newPoint);
        
        omp_set_lock(&treeLock);
        Point *newNode = nodes.insert(Point(newPoint, nn, cost));
        omp_unset_lock(&treeLock);
        
        rewire(newNode);
        return newNode;
    }

    Point* connect(std::vector<double>& point){
        Point* nn = nodes.nearestNeighbor(point);
        do{
            std::vector<double> newPoint = moveStep(nn->coords, point);
            
            if(checkCollision(newPoint)==true)
                return NULL;

            double cost = nn->cost + eulerDistance(nn->coords, newPoint);
            
            omp_set_lock(&treeLock);
            nn = nodes.insert(Point(newPoint, nn, cost));
            omp_unset_lock(&treeLock);
            
            rewire(nn);
        }while(!checkGoalState(nn->coords, point));
        return nn;
    }

    bool checkCollision(std::vector<double>& point){
        return cc.isColliding(point);
    }

    void rewire(Point* newPoint){
        std::vector<Point*> nearPoints = nodes.searchRadius(newPoint->coords, rewireRadius);
        for(auto p: nearPoints){
            double cost = p->cost + eulerDistance(newPoint->coords, p->coords);
            if(cost<newPoint->cost){
                newPoint->parent = p;
                newPoint->cost = cost;
            }
        }
        for(auto p: nearPoints){
            double cost = newPoint->cost + eulerDistance(newPoint->coords, p->coords);
            if(cost<p->cost){
                p->parent = newPoint;
                p->cost = cost;
            }
        }
    }
    
    public:
    
    RRTStar(std::vector<double> startState, std::vector<double> goalCoord, double goalRadius, double stepSize, double rewireRadius, double length):
        nodes(startState.size()),
        cc("map1.dat", length, length)
    {
        nodes.insert(Point(startState));
        this->goalCoord = goalCoord;
        this->goalRadius = goalRadius;
        this->stepSize = stepSize;
        this->rewireRadius = rewireRadius;
        // this->gen = std::mt19937(std::random_device()());
        this->dis = std::uniform_real_distribution<>(-length,length);
        this->dof = startState.size();
    }

    bool doRRTStar(){
        
        long long iter=1;
        Point *final = nodes.getRoot()->getPointAddress();
        bool end_flag = false;
        omp_init_lock(&treeLock);

        std::vector<double> connectPoint(dof,0);  
        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
            {
                while(end_flag==false){

                    Point *temp=NULL;
                    do{
                        temp = extend();
                    }while(temp==NULL);

                    MPI_Send(&(temp->coords[0]), dof, MPI_DOUBLE, (myid+1)%2, 1, MPI_COMM_WORLD);
                }
            }

            #pragma omp section
            {
                MPI_Request req;
                MPI_Irecv(&end_flag, 1, MPI_C_BOOL, (myid+1)%2, 2, MPI_COMM_WORLD, &req);
                while(end_flag==false){     
                        
                    MPI_Recv(&(connectPoint[0]), dof, MPI_DOUBLE, (myid+1)%2, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    final = connect(connectPoint);
                    if(final!=NULL){
                        end_flag = true;
                    }
                    iter++;

                    // MPI_Iallreduce(&end_flag, &end_flag, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
                }
            }
        }
        #pragma omp single
        {
            MPI_Request req;
            MPI_Isend(&end_flag, 1, MPI_C_BOOL, (myid+1)%2, 2, MPI_COMM_WORLD, &req);
            if(final!=NULL)
                MPI_Send(&(connectPoint[0]), dof, MPI_DOUBLE, (myid+1)%2, 3, MPI_COMM_WORLD);
            else{
                std::vector<double> dummy(dof,0);
                MPI_Recv(&(dummy[0]), dof, MPI_DOUBLE, (myid+1)%2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                final = nodes.nearestNeighbor(dummy);
            }
                
            std::ofstream fout("path"+std::to_string(myid)+".dat");
            int stepCount=0;
            for(Point* i=final; i!=NULL; i = i->parent){
                std::cout<<i->coords[0]<<","<<i->coords[1]<<'\n';
                fout<<i->coords[0]<<","<<i->coords[1]<<'\n';
                stepCount++;
            }
            fout.close();

            std::cout<<"Found a path in "<<stepCount<<" steps.\nExpanded "<<nodes.size()<<" nodes.\n";
            nodes.exportTree("nodes"+std::to_string(myid)+".dat");

            omp_destroy_lock(&treeLock);
        }
        return true;
    }
};

std::mt19937 RRTStar::gen = std::mt19937(std::random_device()());

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::vector<double> startState = {-4,-4}, endState = {4,4};
    
    RRTStar r(
        myid==0?startState:endState, 
        myid==0?endState:startState, 
        0.2, 0.2, 0.3, 5
    );

    r.doRRTStar();

    MPI_Finalize();

    return 0;
}
