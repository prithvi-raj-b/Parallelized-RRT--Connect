#include<vector>
#include<random>
#include<cmath>
#include<iostream>
#include<fstream>
#include<mpi.h>

#include "collision.cpp"

#define N 2
#define goal_tol 0.1 //tolerance
#define default_node_step_size 0.1 //goal tol taken to be equal to this, or lesser!!!
#define node_radius 0.3 //radius of rewiring for rrt*
#define id 1

//RRT_Star, with obstacle checking

CollisionChecker coll("map1.dat",5,5);

class node{
    public:
    
    std::vector<double> point;
    double cost;
    long parent;  //every node has just one parent! this gets changed everytime we do the RRT star step  

    node(std::vector<double> point): point(point), cost(0), parent(-1) {
    }

    node(std::vector<double> point, double cost, long parent): point(point), cost(cost), parent(parent) {
    }

    node(): point(N,0), cost(0), parent(-1){
    }

    void update(node& new_node)
    {
        point = new_node.point;
        cost = new_node.cost;
        parent = new_node.parent;
    }
};

double norm(std::vector<double>& point1, std::vector<double>& point2)
{
    double result = 0.0;
    long unsigned int i;

    for(i=0; i<point1.size(); i++)
    {
        result += (point1.at(i) - point2.at(i))*(point1.at(i) - point2.at(i));
    }
    result = sqrt(result);
    return result;
}

// create random point in main program.
// function1 to compare distance to all values in graph, return graph index of potential parent to spawn link from (only potential because we need to collision check). 
long candidate_link(std::vector<node>& graph, std::vector<double>& random_point)
{
    double min_distance;
    long candidate_link_index = 0;

    long unsigned int i;

    min_distance = norm(graph.at(0).point,random_point);

    //std::cout<<"distance : "<<min_distance<<std::endl;

    //loop to find actual values

    for(i=1; i<graph.size(); i++)
    {
        //std::cout<<norm(graph.at(i).point,random_point)<<std::endl;
        if(min_distance >= norm(graph.at(i).point,random_point)) 
        {
            min_distance = norm(graph.at(i).point,random_point);
            //std::cout<<"min_distance = "<<min_distance<<std::endl;
            candidate_link_index = i;
        }
    }

    //std::cout<<"candidate_link_index = "<<candidate_link_index<<std::endl;
    return candidate_link_index;
}

//-------------------------
//function to spawn a link between random_point and chosen graph point indicated by link_index. It spawns only step_distance if 
long long spawn_node(long link_index, std::vector<node>& graph, std::vector<double>& random_point)
{
    long unsigned int i;

    double step_distance = default_node_step_size; 
    
    std::vector<double> slope(N);

    double dist = norm(graph.at(link_index).point,random_point);
    //std::cout<<"distance :"<<dist<<std::endl;

    for(i=0; i<random_point.size();i++)
    {
        slope[i] = (random_point.at(i)-graph.at(link_index).point.at(i))/dist; //std::cout<<"ok!"<<std::endl;
    }

    node temp(graph.at(link_index));

    if(dist < step_distance) 
        step_distance = dist;

    for(i=0; i<random_point.size();i++) 
        temp.point[i] += slope.at(i)*step_distance;
    
    if(coll.isColliding(temp.point)) 
        return -1; //if collision, don't spawn node
    
    temp.cost += step_distance;

    temp.parent = link_index; //set by default. Can be updated if cheaper connect found!

    long long res;
    #pragma omp critical
    {
        graph.push_back(temp);  //this adjoins a node along the direction picked by the random_point along a distance of only 0.1
        res = graph.size()-1; //return the index of the new node spawned
    }
    return res;
}

//------------------
// function2 to take in previous random point and candidate parent to collision check. If returns true(no collision) update parent and cost of new node, and push node into the graph, then proceed. Else go back and choose a different random point


//--------------------------------------------------
// STAR STEP : Now check some proximity of new node for cheaper spawns, i.e. a candidate spawn for which cost of new_node is lesser than the present value
// Do collision checking for the line between the new node and any point within proximity. Then check if this point makes the new node cheaper than it presently is. If yes, update the index to be returned and keep checking for cheaper spawn points

//this returns the graph (vector) index of the cheapest spawn for the given random node
long cheapest_spawn_index(long node_spawn_index, std::vector<double>& new_node_point, std::vector<node>& graph, std::vector<long>& addresses_in_node_radius) 
{
    double dist,new_cost;
    long unsigned int i;

    node cheapest_spawn(graph.at(node_spawn_index));

    long cheapest_connect_index = node_spawn_index;

    for(i=0; i<graph.size(); i++)
    {
        dist = norm(graph[i].point,new_node_point);
        new_cost = dist + graph[i].cost;

        if(dist < node_radius)    
        {
            addresses_in_node_radius.push_back(i); //store all addresses in node radius into this vector
            
            if(new_cost < cheapest_spawn.cost + dist  && !coll.isColliding(new_node_point)) 
            {
                cheapest_spawn.update(graph[i]);
                cheapest_connect_index = i;
            }
        }
    }
    //std::cout<<"cheapest_spawn_index = "<<cheapest_spawn_index<<std::endl;
    return cheapest_connect_index;
}

int main (int argc, char **argv)
{
    std::vector<node> graph(1);

    node final_node;

    node start({-4,-4},0,-1); 
    node finish({4,4},0,-1); 

    std::vector<double> end(2);

    long cnt=0;
    std::vector<double> newest_of_other(2);

    std::uniform_real_distribution<double> unif(-4.999,4.999);  //5 x 5 grid say
    std::default_random_engine re = std::default_random_engine(std::random_device()());

    int  myid; //nprocs = 2 fixed. 0 grows from the start, 1 grows from finish

    int end_flag;

    

    int tag1=1;
    long long spawn_int = -1;
    if(myid == id) std::cout<<"Reached0: "<<cnt<<std::endl;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    //save the end goal for each graph
    if(myid == 0) 
    {
        graph[0].update(start);
        end.at(0) = finish.point.at(0); end.at(1) = finish.point.at(1); //hardcode N=2
    }  
    else 
    {
        graph[0].update(finish);
        end.at(0) = start.point.at(0); end.at(1) = start.point.at(1);  //hardcode N=2
    }

    end_flag=0;
    if(myid == id) std::cout<<"Reached0: "<<cnt<<std::endl;

    #pragma omp parallel sections num_threads(2) private(re) shared(std::cout,graph, end_flag, end, cnt,unif,tag1,myid,final_node)
    {
        #pragma omp section 
        {
            //pick a random point and spawn a node on the closest graph member to it, in the direction towards this random point.
            long long spawn_int;
            bool start_pop;
            std::vector<double> random_point(2);
            std::vector<long> addresses_in_node_radius;
            long link_index=0;
            while(cnt <= 1000 && !end_flag)
            {
                do
                {
                    random_point.at(0) = unif(re);
                    random_point.at(1) = unif(re);

                    link_index = candidate_link(graph,random_point); //returns graph address of nearest link
                    spawn_int = spawn_node(link_index,graph,random_point);

                    if(myid==id) std::cout<<"hanging? "<<cnt<<std::endl;

                }while( spawn_int==-1 );  //exit loop only if no collision

                if(myid==id) std::cout<<"Reached1 "<<cnt<<std::endl;
            
                addresses_in_node_radius.clear();

                graph[spawn_int].parent =  cheapest_spawn_index(link_index,graph[spawn_int].point,graph,addresses_in_node_radius); 

                if(myid==id) std::cout<<"Reached2 "<<cnt<<std::endl;
                addresses_in_node_radius.push_back(spawn_int); //last element added is the address of the new node! skip this address in rewiring!
                for(long unsigned int i=0; i<addresses_in_node_radius.size()-1; i++) 
                {
                    if(graph[addresses_in_node_radius[i]].cost > norm(graph[addresses_in_node_radius[i]].point,graph[spawn_int].point) + graph[spawn_int].cost) 
                        graph[addresses_in_node_radius[i]].parent = spawn_int;
                
                } //rewiring possible TO NEWEST NODE ONLY, if cheaper.

                //now we have grown one RRT* node in both graphs. To connect nearest neighbour of OTHER graph to NEW NEIGHBOUR of this graph.

                MPI_Send(&graph[spawn_int].point[0],2,MPI_DOUBLE,(myid+1)%2,tag1,MPI_COMM_WORLD); //send most recent member of this graph to other process
                MPI_Barrier(MPI_COMM_WORLD);
                if(myid==id) std::cout<<"Reached3 "<<cnt<<std::endl;
                cnt++;
            }
        }

        #pragma omp section
        {
            long link_index=0;
            while(cnt <= 1000 && !end_flag)
            {
                MPI_Recv(&newest_of_other[0],2,MPI_DOUBLE,(myid+1)%2,tag1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); //receive the above point in a vector called "newest_of_other"

                MPI_Barrier(MPI_COMM_WORLD); //to make sure send-recieves are synced
                
                link_index = candidate_link(graph,newest_of_other); // finds address of node in this process's graph, which is closest to newest_of_other

                
                //the following block of code grow the graphs until an obstacle is reached or they get close enough to a point.

                if(myid == id) std::cout<<"Reached0: "<<cnt<<std::endl;
                spawn_int = spawn_node(link_index,graph,newest_of_other);
                if(myid == id) std::cout<<"Reached1: "<<spawn_int<<std::endl;
                while(spawn_int != -1)
                {    
                    if(myid == id) std::cout<<"am i infinite : "<<cnt<<std::endl;
                    
                    link_index = spawn_int;
                    spawn_int = spawn_node(link_index,graph,newest_of_other);  

                    if(norm(newest_of_other,graph[link_index].point)< goal_tol || norm(end,graph[link_index].point) < goal_tol)
                    {
                        final_node.parent = link_index;

                        if(norm(newest_of_other,graph[link_index].point)< goal_tol) 
                        {
                            final_node.point[0] = newest_of_other[0];
                            final_node.point[1] = newest_of_other[1];
                            final_node.cost = graph[link_index].cost + norm(newest_of_other,graph[link_index].point);
                        }
                        
                        else
                        {
                            final_node.point[0] = end[0];
                            final_node.point[1] = end[1];
                            final_node.cost = graph[link_index].cost + norm(end,graph[link_index].point);              
                        }
                        end_flag=1;
                        if(myid == id) std::cout<<"Completed!"<<std::endl;
                        break;
                    }
                }

                if(end_flag) graph.push_back(final_node);

                //send end_flag to both procs. 
                MPI_Allreduce(&end_flag,&end_flag,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);
                cnt++;
            }
        }
    } 

        if(myid == id) std::cout<<"not done "<<cnt<<std::endl;

    // }
    //there is a fundamental difference between this approach and the RRT connect approach given by Kuffner and LaValle, 
    //they do this serially by SWAPPING the graphs, whereas we identify the closest points before growing parallely and grow then simultaneously.

    //write stuff into files for plotting

    #pragma omp single
    if(myid == 0)
    {
        std::ofstream myfile; 
        myfile.open("graph_start.txt");
        for (unsigned long i = 1; i < graph.size(); i++)
        {
            myfile<<graph[i].point[0]<<" "<<graph[i].point[1]<<" "<<graph[graph[i].parent].point[0]<<" "<<graph[graph[i].parent].point[1]<<std::endl;
        }
        myfile.close();
        
    }
    else 
    {
        std::ofstream myfile; 
        myfile.open("graph_finish.txt");
        for (unsigned long i = 1; i < graph.size(); i++)
        {
            myfile<<graph[i].point[0]<<" "<<graph[i].point[1]<<" "<<graph[graph[i].parent].point[0]<<" "<<graph[graph[i].parent].point[1]<<std::endl;
        }
        myfile.close();
    }

    #pragma omp single
    MPI_Finalize();
    
}


