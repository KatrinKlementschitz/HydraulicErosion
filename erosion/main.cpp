#include <SFML/Graphics.hpp>
#include "PerlinNoise.h"
#include <iostream>
#include <random>

#define cellsize 2
#define len 400
#define inertia 0.45
#define M_PI 3.14159265358979323846

double initialWaterVolume = 1;
double initialSpeed = 1;
int maxDropletLifetime = 120;
double erodeSpeed = .9f;
double depositSpeed = 1.f;
double sedimentCapacityFactor = 4;
double evaporateSpeed = .05f;
double minSlope = 0.0001;
double gravity = 9.81;
double minSedimentCapacity = .01f;
double stepFac = 1.f;

template <class T>
T max(T a, T b) {
    return (a > b ? a : b);
}

template <class T>
T min(T a, T b) {
    return (a < b ? a : b);
}

struct gradient
{
    double gradientX;
    double gradientY;
    double Height;
};

struct dropletMove
{
    double x;
    double y;
    double gradientx;
    double gradienty;
    double length;
};

struct dropletPoint
{
    double x;
    double y;
    double gradientx;
    double gradienty;
    double length;
    double deltaHeight;
    double newHeight;
};

class Grid
{
    std::vector<double> arr;
    int rows;
    int columns;
public:
    Grid(int rows, int columns)
    {
        arr = std::vector<double>(rows*columns,0);
        this->rows = rows;
        this->columns = columns;
    }

    double GetVal(int row, int column)
    {
        if(row < rows && column < columns)
			return arr[column * rows + row];
    }

    void SetVal(int row, int column, double val)
    {
        if (row < rows && column < columns)
	        if(!isnan(val) && !isinf(val))
				arr[column * rows + row] = val;
    }

    void AddVal(int row, int column, double val)
    {
        /*if(arr[column * rows + row] + val < 0)
            return;*/
        if (row < rows && column < columns)
	        if (!isnan(val) && !isinf(val))
	            arr[column * rows + row] += val;
    }

    void MultVal(int row, int column, double val)
    {
        if (row < rows && column < columns)
	        if (!isnan(val) && !isinf(val))
	            arr[column * rows + row] *= val;
    }
};

Grid grid(len/cellsize, len/cellsize);

int iter = 0;

void NormalizeGrid()
{
    double min = 100;
    double max = -100;

    for (int i = 0; i < len / cellsize; i++)
    {
        for (int j = 0; j < len / cellsize; j++)
        {
            if(!isnan(grid.GetVal(i,j)))
	            if (grid.GetVal(i, j) > max)
	                max = grid.GetVal(i, j);
	            if (grid.GetVal(i, j) < min)
	                min = grid.GetVal(i, j);
        }
    }

    double mult = (max - min);
    for (int i = 0; i < len / cellsize; i++)
    {
        for (int j = 0; j < len / cellsize; j++)
        {
            double temp = grid.GetVal(i, j);
            double temp1 = temp - min;
            double temp2 = temp1 / mult;

            grid.SetVal(i, j, temp1);
            if (isnan(temp))
                grid.SetVal(i, j, 0);
            grid.SetVal(i, j, temp2);
        }
    }
}

void Blur(double f1)
{
    double f2 = (3.0 - f1) / 2.0;
    for (int j = 1; j < len / cellsize - 1; j++)
    {
        for (int i = 1; i < len / cellsize - 1; i++)
        {
            double k = (grid.GetVal(i - 1, j)*f2 + f1*grid.GetVal(i, j) + f2*grid.GetVal(i + 1, j)) / 3.0;
            grid.SetVal(i, j, k);
        }
        for (int i = 1; i < len / cellsize - 1; i++)
        {
            double k = (grid.GetVal(j, i - 1)*f2 + f1*grid.GetVal(j, i) + f2*grid.GetVal(j, i + 1)) / 3.0;
            grid.SetVal(j, i, k);
        }
    }
}

void Brush(int centreX, int centreY, int radius, double val)
{
    double weightSum = 0;
    int addIndex = 0;
    for (int y = -radius; y <= radius; y++) {
        for (int x = -radius; x <= radius; x++) {
	        const float sqrDst = x * x + y * y;
            if (sqrDst < radius * radius) {
	            const int coordX = centreX + x;
	            const int coordY = centreY + y;

	            const float weight = 1 - sqrt(sqrDst) / radius;
                weightSum += weight;
                addIndex++;

                if (coordX < len / cellsize && coordX >= 0 && coordY < len / cellsize && coordY >= 0)
                    grid.AddVal(coordX, coordY, - weight * val);
            }
        }
    }
}

double calcCarryCapacity(double heightDiff, double minSlope, double velocity, double water, double capacity)
{
    return max(abs(heightDiff), minSlope) * velocity * water * capacity;
}

double calcHeightDiff(double HeightOld, double HeightNew)
{
    return HeightNew - HeightOld;
}

double calcVelocity(double velocity, double HeightDiff, double gravity)
{
    return  sqrt(velocity * velocity + HeightDiff * gravity);
}

double calcWater(double water, double evaporation)
{
    return  water * (1 - evaporation);
}

gradient calcGradient(double posX, double posY)
{
    gradient g;
    const int coordX = static_cast<int>(posX);
    const int coordY = static_cast<int>(posY);

    double x = posX - coordX;
    double y = posY - coordY;

    double heightNW, heightNE, heightSW, heightSE;

    heightNW = grid.GetVal(coordX - 1, coordY - 1);
    heightNE = grid.GetVal(coordX + 1, coordY - 1);
    heightSW = grid.GetVal(coordX - 1, coordY + 1);
    heightSE = grid.GetVal(coordX + 1, coordY + 1);

    const double gradientX = (heightNE - heightNW) * (1 - y) + (heightSE - heightSW) * y;
    const double gradientY = (heightSW - heightNW) * (1 - x) + (heightSE - heightNE) * x;

    const double height = heightNW * (1 - x) * (1 - y) + heightNE * x * (1 - y) + heightSW * (1 - x) * y + heightSE * x * y;

    g.gradientX = gradientX;
    g.gradientY = gradientY;
    g.Height = height;
    
    return g;
}

dropletMove CalcDropletMove(gradient g)
{
    dropletMove p;
    double dirX = 0;
    double dirY = 0;
    dirX = dirX * inertia - g.gradientX * (1 - inertia);
    dirY = dirY * inertia - g.gradientY * (1 - inertia);
    const double length = sqrt(dirX * dirX + dirY * dirY);
    dirX = dirX / length / stepFac;
    dirY = dirY / length / stepFac;
    p.gradientx = g.gradientX;
    p.gradienty = g.gradientY;
    p.x = dirX;
    p.y = dirY;
    p.length = length;
    return p;
}

dropletPoint CalcNewDropletPoint(double cx, double cy)
{
    dropletPoint pr;
    gradient g = calcGradient(cx, cy);
    dropletMove p = CalcDropletMove(g);

    const int gridx = cx;
    const int gridy = cy;

    const double cellOffsetX = cx - gridx;
    const double cellOffsetY = cy - gridy;
    cx = p.x + cx;
    cy = p.y + cy;

    const double newHeight = calcGradient(cx, cy).Height;

    const double deltaHeight = calcHeightDiff(g.Height, newHeight);

    pr.x = cx;
    pr.y = cy;
    pr.length = p.length;
    pr.deltaHeight = deltaHeight;
    pr.newHeight = newHeight;
    pr.gradientx = p.gradientx;
    pr.gradienty = p.gradienty;
    return  pr;
}

void PlaceSediment(double amountToDeposit, int gridx, int gridy, double cellOffsetX, double cellOffsetY, double min, double max, double current)
{
    grid.AddVal(gridx, gridy, amountToDeposit * (1 - cellOffsetX) * (1 - cellOffsetY));// * (min / current));
    grid.AddVal(gridx + 1, gridy, amountToDeposit * cellOffsetX * (1 - cellOffsetY));// *(min / grid.GetVal(gridx + 1, gridy)));
    grid.AddVal(gridx, gridy + 1, amountToDeposit * (1 - cellOffsetX) * cellOffsetY);// * (min / grid.GetVal(gridx, gridy + 1)));
    grid.AddVal(gridx + 1, gridy + 1, amountToDeposit * cellOffsetX * cellOffsetY);// * (min / grid.GetVal(gridx + 1, gridy + 1)));
}

void Erode(double mx, double my) {
    erodeSpeed = .01f;
    double cx = mx / cellsize;
    double cy = my / cellsize;
    int i = 0;

    double speed = initialSpeed;
    double water = initialWaterVolume;
    double sediment = 0;

    double oldgx = 0;
    double oldgy = 0;

    while (i < maxDropletLifetime)
    {
        const int cx_int = static_cast<int>(cx);
        const int cy_int = static_cast<int>(cy);

        if (cx_int > 0 && cx_int < len / cellsize - 1 && cy_int > 0 && cy_int < len / cellsize - 1 && !isnan(cx) && !isnan(cy))
        {
            i++;
            const double cellOffsetX = cx - cx_int;
            const double cellOffsetY = cy - cy_int;
            dropletPoint p = CalcNewDropletPoint(cx, cy);
            cx = p.x;
            cy = p.y;

            if (isnan(cx) || isnan(cy)){
                return;
            }
            const double deltaHeight = p.deltaHeight;

            if (abs(deltaHeight) > 10) {
                return;
            }

            if (isnan(speed))
                speed = 0.0001;


            //const double sedimentCapacity = max(-deltaHeight * speed * water * sedimentCapacityFactor, minSedimentCapacity);
            //const double sedimentCapacity = max(-deltaHeight * speed * water * sedimentCapacityFactor, minSedimentCapacity);

            const double sedimentCapacity = calcCarryCapacity(deltaHeight, minSlope, speed, water, sedimentCapacityFactor);

            const double loose_sediment = ((oldgx - p.gradientx) + (oldgy - p.gradienty)) / (p.gradientx+p.gradienty);

            oldgx = p.gradientx;
            oldgy = p.gradienty;


            if (!isnan(cx) && !isnan(cy))
            {
                const double minimum = min(min(grid.GetVal(cx_int, cy_int), grid.GetVal(cx_int + 1, cy_int)), min(grid.GetVal(cx_int, cy_int + 1), grid.GetVal(cx_int + 1, cy_int + 1)));
                const double maximum = max(max(grid.GetVal(cx_int, cy_int), grid.GetVal(cx_int + 1, cy_int)), max(grid.GetVal(cx_int, cy_int + 1), grid.GetVal(cx_int + 1, cy_int + 1)));
                const double current = grid.GetVal(cx_int, cy_int);

                if (loose_sediment > 1.0)
                {
                    erodeSpeed /= 2.0;
                    //PlaceSediment(abs(sediment), cx_int, cy_int, cellOffsetX, cellOffsetY, minimum, maximum, current);
                    return;
                }

                if (p.length < 0.002 || i == maxDropletLifetime-1)
                {
                    Brush(cx_int, cy_int, 2, -abs(sediment));
                	return;
                }

                if (sediment > sedimentCapacity || deltaHeight > 0) {
                    // If moving uphill (deltaHeight > 0) try fill up to the current height, otherwise deposit a fraction of the excess sediment
                    double amountToDeposit = (deltaHeight > 0)
                        ? min(deltaHeight, sediment)
                        : (sediment - sedimentCapacity) * depositSpeed;
                    amountToDeposit = abs(amountToDeposit);
                    sediment -= amountToDeposit;

                    PlaceSediment(amountToDeposit, cx_int, cy_int, cellOffsetX, cellOffsetY, minimum, maximum, current);
                }
                else
                {
                    double amountToErode = min((sedimentCapacity - sediment) * erodeSpeed, -deltaHeight);
                    amountToErode = abs(amountToErode);
                    sediment += amountToErode;

					if (!isnan(amountToErode)) {
                        if (amountToErode < p.newHeight){
                            Brush(cx_int, cy_int, 2, amountToErode);
                        }
                    }
                }

                speed = calcVelocity(speed, deltaHeight, gravity);
                water = calcWater(water, evaporateSpeed);
            }
            if (p.length < 0.002) { return; }
        }
        else { return; }
    }
}

//Red deposit, Blue erode
void DrawPath(sf::RenderWindow &window, double mx, double my)
{
    sf::Vertex line[] =
    {
        sf::Vertex(sf::Vector2f(10.f, 10.f)),
        sf::Vertex(sf::Vector2f(150.f, 150.f))
    };
    line[0].color.r = 255;
    line[0].color.g = 0;
    line[0].color.b = 0;
    line[1].color.r = 255;
    line[1].color.g = 0;
    line[1].color.b = 0;

    double cx = mx / cellsize;
    double cy = my / cellsize;
    int i = 0;

    double speed = initialSpeed;
    double water = initialWaterVolume;
    double sediment = 0;
    int oldgx = 0;
    int oldgy = 0;
    int samecell = 0;

    while (i<maxDropletLifetime)
    {
        if ((int)cx > 0 && (int)cx < len / cellsize - 1 && (int)cy > 0 && (int)cy < len / cellsize - 1)
        {
            

            i++;
            gradient g = calcGradient(cx, cy);


            double dirX = 0;
            double dirY = 0;
            dirX = (dirX * inertia - g.gradientX * (1 - inertia));
            dirY = (dirY * inertia - g.gradientY * (1 - inertia));
            double length = sqrt(dirX * dirX + dirY * dirY);
            dirX = dirX / length / stepFac;
            dirY = dirY / length / stepFac;
            line[0].position.x = cx * cellsize;
            line[0].position.y = cy * cellsize;
            cx = dirX + cx;
            cy = dirY + cy;
            line[1].position.x = cx*cellsize;
            line[1].position.y = cy*cellsize;

            if(isnan(cx) || isnan(cy))
                break;

            double newHeight = calcGradient(cx, cy).Height;
            double deltaHeight = calcHeightDiff(g.Height, newHeight);
            

            if (isnan(speed))
                speed = 0.1;

            const double sedimentCapacity = calcCarryCapacity(deltaHeight, minSlope, speed, water, sedimentCapacityFactor);

            double cellOffsetX = cx - (int)cx;
            double cellOffsetY = cy - (int)cy;
            if (!isnan(cx) && !isnan(cy))
            {

                if (sediment > sedimentCapacity || deltaHeight > 0) {
                    // If moving uphill (deltaHeight > 0) try fill up to the current height, otherwise deposit a fraction of the excess sediment
                    double amountToDeposit = (deltaHeight > 0) ? min(deltaHeight, sediment) : (sediment - sedimentCapacity) * depositSpeed;
                    sediment -= amountToDeposit;
                    line[0].color.r = 255;
                    line[1].color.r = 255;
                    line[0].color.b = 0;
                    line[1].color.b = 0;
                }
                else
                {
                    float amountToErode = min((sedimentCapacity - sediment) * erodeSpeed, -deltaHeight);
                    sediment += amountToErode;
                    line[0].color.b = 255;
                    line[1].color.b = 255;
                    line[0].color.r = 0;
                    line[1].color.r = 0;
                }
            }

            speed = calcVelocity(speed, deltaHeight, gravity);
            water = calcWater(water, evaporateSpeed);

            window.draw(line, 2, sf::Lines);
            if(length<0.002)
                break;
        }
        else { break; }
    }
}

int main()
{
    sf::RenderWindow window(sf::VideoMode(len, len), "Hydraulic Erosion!");
    PerlinNoise pn;
    sf::Vertex line[] =
    {
        sf::Vertex(sf::Vector2f(10.f, 10.f)),
        sf::Vertex(sf::Vector2f(150.f, 150.f))
    };
    line[0].color.r = 255;
    line[0].color.g = 0;
    line[0].color.b = 0;
    line[1].color.r = 0;
    line[1].color.g = 0;
    line[1].color.b = 255;

    srand(time(NULL));
    int r = rand() % 1000;
    window.setMouseCursorVisible(false);
    window.setVerticalSyncEnabled(false);

    sf::RectangleShape rectangle(sf::Vector2f(cellsize, cellsize));
    rectangle.setFillColor(sf::Color::Blue);
    rectangle.setPosition(0, 0);
    rectangle.setOutlineThickness(0);

    sf::RectangleShape point(sf::Vector2f(1, 1));
    point.setFillColor(sf::Color::Blue);
    point.setPosition(0, 0);
    point.setOutlineThickness(0);

    double fac = 0.01;

    for (int i = 0; i < len/cellsize; i++)
    {
        for (int j = 0; j < len/cellsize; j++)
        {
            double n = pn.noise(i * fac, j * fac, r);
            n += pn.noise(i * fac * 2, j * fac * 2, r + 10) * 0.5;
            n += pn.noise(i * fac * 4, j * fac * 4, r + 20) * 0.25;
            n += pn.noise(i * fac * 8, j * fac * 8, r + 30) * 0.125;
            grid.SetVal(i, j, n);
        }
    }
    NormalizeGrid();
    double lower_bound = 0;
    double upper_bound = len;
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine re;

    for (int i = 0; i < 100000; i++)
    {
        double rx = unif(re);
        double ry = unif(re);

        Erode(rx, ry);
    }
    //NormalizeGrid();
	//Blur(2.7);

    double maxX = 0;
    double maxY = 0;
    double maxH = 0;

    for (int i = 1; i < len / cellsize; i++)
    {
        for (int j = 1; j < len / cellsize; j++)
        {
            gradient x = calcGradient(i, j);
            if (abs(x.gradientX) > maxX && abs(x.gradientX) < 10)
                maxX = abs(x.gradientX);
        	if (abs(x.gradientY) > maxY && abs(x.gradientY) < 10)
                maxY = abs(x.gradientY);
        	if (abs(x.Height) > maxH && abs(x.Height))
                maxH = x.Height;
        }
    }


    sf::Color a;
    a.b = 127;
    a.r = 127;
    a.g = 127;

    int mx = 0;
    int my = 0;
    iter = 0;

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();

            if(event.type == sf::Event::MouseMoved){
                mx = event.mouseMove.x;
                my = event.mouseMove.y;
            }

            if (event.type == sf::Event::MouseButtonPressed) {
                if(event.mouseButton.button == sf::Mouse::Right){
	                for (int i = 0; i < 7000; i++)
	                {
	                    double rx = unif(re);
	                    double ry = unif(re);

	                    Erode(rx, ry);
                        NormalizeGrid();
	                }
                }
                else
                {
                    std::cout << grid.GetVal(mx/cellsize, my/cellsize) << std::endl;
                }
            }
        }

        window.clear();

        for (int i = 0; i < len/cellsize; i++)
        {
            for (int j = 0; j < len/cellsize; j++)
            {
                double n = grid.GetVal(i, j);
                gradient g = calcGradient(i, j);
                a.b = g.Height / maxH * 255.0;
                a.r = a.b;
                a.g = a.r;
                rectangle.setFillColor(a);
                rectangle.setPosition(i * cellsize, j * cellsize);
                window.draw(rectangle);
            }
        }

        double rx = unif(re);
        double ry = unif(re);
        DrawPath(window, mx, my);
        //window.setTitle(std::to_string(iter));
        //for (int i = 0; i < 100; i++)
        //{
        //    double rx = unif(re);
        //    double ry = unif(re);

        //    //DrawPath(window, rx, ry);
        //    Erode(rx, ry);
        //}
        //iter++;
        window.display();
    }

    return 0;
}