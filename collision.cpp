#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<cmath>

class CollisionChecker{
    std::string fname;
    std::vector<std::vector<char>> map;
    double length, width;

    static char atoi(char c){
        if(c=='0')
            return 0;
        else if(c=='1')
            return 1;
        return 1;
    }

    public:

    CollisionChecker(std::string fname, double length, double width): fname(fname), length(length), width(width) {

        std::ifstream fin(fname);
        int n,m;
        fin>>n>>m;
        map.resize(n, std::vector<char>(m));
        for(auto &row: map)
            for(auto &cell: row){
                fin>>cell;
                cell = atoi(cell);
            }
        fin.close();
    }
    
    bool isColliding(std::vector<double> &point){
        // TODO: interpolated checking
        int y = round((point[0]+length)/(2*length)*map.size());
        int x = round((point[1]+width)/(2*width)*map[0].size());
        if(x<0 or y<0 or x>=map[0].size() or y>=map.size())
            return true;
        return map[x][y]==1;
    }
};

// int main(){
//     CollisionChecker cc("map1.dat", 5, 5);
//     std::vector<std::vector<double>> point = {{-1.8,-4},{-3,-4},{-1.8,3},{-6,-6}};
//     for(auto p: point)
//         std::cout<<cc.isColliding(p)<<std::endl;
//     return 0;
// }
