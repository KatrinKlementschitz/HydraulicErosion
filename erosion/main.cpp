#include <SFML/Graphics.hpp>
#include "PerlinNoise.h"
#include <iostream>
#include <random>
#include  "erosion.h"

#define cellsize 2
#define len 400
#define M_PI 3.14159265358979323846


//void Blur(double f1)
//{
//    double f2 = (3.0 - f1) / 2.0;
//    for (int j = 1; j < len / cellsize - 1; j++)
//    {
//        for (int i = 1; i < len / cellsize - 1; i++)
//        {
//            double k = (grid.GetVal(i - 1, j)*f2 + f1*grid.GetVal(i, j) + f2*grid.GetVal(i + 1, j)) / 3.0;
//            grid.SetVal(i, j, k);
//        }
//        for (int i = 1; i < len / cellsize - 1; i++)
//        {
//            double k = (grid.GetVal(j, i - 1)*f2 + f1*grid.GetVal(j, i) + f2*grid.GetVal(j, i + 1)) / 3.0;
//            grid.SetVal(j, i, k);
//        }
//    }
//}


//Red deposit, Blue erode
//void DrawPath(sf::RenderWindow &window, double mx, double my)
//{
//    sf::Vertex line[] =
//    {
//        sf::Vertex(sf::Vector2f(10.f, 10.f)),
//        sf::Vertex(sf::Vector2f(150.f, 150.f))
//    };
//    line[0].color.r = 255;
//    line[0].color.g = 0;
//    line[0].color.b = 0;
//    line[1].color.r = 255;
//    line[1].color.g = 0;
//    line[1].color.b = 0;
//
//    double cx = mx / cellsize;
//    double cy = my / cellsize;
//    int i = 0;
//
//    double speed = initialSpeed;
//    double water = initialWaterVolume;
//    double sediment = 0;
//    int oldgx = 0;
//    int oldgy = 0;
//    int samecell = 0;
//
//    while (i<maxDropletLifetime)
//    {
//        if ((int)cx > 0 && (int)cx < len / cellsize - 1 && (int)cy > 0 && (int)cy < len / cellsize - 1)
//        {
//            
//
//            i++;
//            gradient g = calcGradient(cx, cy);
//
//
//            double dirX = 0;
//            double dirY = 0;
//            dirX = (dirX * inertia - g.gradientX * (1 - inertia));
//            dirY = (dirY * inertia - g.gradientY * (1 - inertia));
//            double length = sqrt(dirX * dirX + dirY * dirY);
//            dirX = dirX / length / stepFac;
//            dirY = dirY / length / stepFac;
//            line[0].position.x = cx * cellsize;
//            line[0].position.y = cy * cellsize;
//            cx = dirX + cx;
//            cy = dirY + cy;
//            line[1].position.x = cx*cellsize;
//            line[1].position.y = cy*cellsize;
//
//            if(isnan(cx) || isnan(cy))
//                break;
//
//            double newHeight = calcGradient(cx, cy).Height;
//            double deltaHeight = calcHeightDiff(g.Height, newHeight);
//            
//
//            if (isnan(speed))
//                speed = 0.1;
//
//            const double sedimentCapacity = calcCarryCapacity(deltaHeight, minSlope, speed, water, sedimentCapacityFactor);
//
//            double cellOffsetX = cx - (int)cx;
//            double cellOffsetY = cy - (int)cy;
//            if (!isnan(cx) && !isnan(cy))
//            {
//
//                if (sediment > sedimentCapacity || deltaHeight > 0) {
//                    // If moving uphill (deltaHeight > 0) try fill up to the current height, otherwise deposit a fraction of the excess sediment
//                    double amountToDeposit = (deltaHeight > 0) ? min(deltaHeight, sediment) : (sediment - sedimentCapacity) * depositSpeed;
//                    sediment -= amountToDeposit;
//                    line[0].color.r = 255;
//                    line[1].color.r = 255;
//                    line[0].color.b = 0;
//                    line[1].color.b = 0;
//                }
//                else
//                {
//                    float amountToErode = min((sedimentCapacity - sediment) * erodeSpeed, -deltaHeight);
//                    sediment += amountToErode;
//                    line[0].color.b = 255;
//                    line[1].color.b = 255;
//                    line[0].color.r = 0;
//                    line[1].color.r = 0;
//                }
//            }
//
//            speed = calcVelocity(speed, deltaHeight, gravity);
//            water = calcWater(water, evaporateSpeed);
//
//            window.draw(line, 2, sf::Lines);
//            if(length<0.002)
//                break;
//        }
//        else { break; }
//    }
//}

int main()
{
    sf::RenderWindow window(sf::VideoMode(len, len), "Hydraulic Erosion!");
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

    //window.setMouseCursorVisible(false);
    window.setVerticalSyncEnabled(false);

    sf::RectangleShape rectangle(sf::Vector2f(cellsize, cellsize));
    rectangle.setFillColor(sf::Color::Blue);
    rectangle.setPosition(0, 0);
    rectangle.setOutlineThickness(0);

    sf::RectangleShape point(sf::Vector2f(1, 1));
    point.setFillColor(sf::Color::Blue);
    point.setPosition(0, 0);
    point.setOutlineThickness(0);

    double maxH = 0;

    Erosion::Config config;
    Erosion erosion(len/cellsize, len/cellsize, 0, config);

    for (int i = 1; i < len / cellsize; i++)
    {
        for (int j = 1; j < len / cellsize; j++)
        {
            double val = erosion.grid->GetVal(i, j);
        	if (abs(val) > maxH && abs(val))
                maxH = val;
        }
    }


    sf::Color a;
    a.b = 127;
    a.r = 127;
    a.g = 127;

    int mx = 0;
    int my = 0;
    int iter = 0;

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
                    erosion.Erode(10000);
                }
                else
                {
                    std::cout << erosion.grid->GetVal(mx/cellsize, my/cellsize) << std::endl;
                }
            }
        }

        window.clear();

        for (int i = 0; i < len/cellsize; i++)
        {
            for (int j = 0; j < len/cellsize; j++)
            {
                double n = erosion.grid->GetVal(i, j);
                a.b = n / maxH * 255.0;
                a.r = a.b;
                a.g = a.r;
                rectangle.setFillColor(a);
                rectangle.setPosition(i * cellsize, j * cellsize);
                window.draw(rectangle);
            }
        }

        erosion.DrawPath(window, mx, my, cellsize);
        erosion.Erode(1000);
        //window.setTitle(std::to_string(iter));
        //iter++;
        window.display();
    }

    return 0;
}