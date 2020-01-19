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
