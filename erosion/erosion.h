#pragma once
#include <vector>
#include <random>
#include <SFML/Graphics/RenderWindow.hpp>

class Grid
{
    std::vector<double> arr;
    int rows;
    int columns;
    int length;
public:
    Grid(int rows, int columns)
    {
        arr = std::vector<double>(rows * columns, 0);
        this->rows = rows;
        this->columns = columns;
        length = rows * columns;
    }

    double GetVal(int row, int column)
    {
        if (row < rows && column < columns)
            return arr[column * rows + row];
    }

    void SetVal(int row, int column, double val)
    {
        if (row < rows && column < columns)
            if (!isnan(val) && !isinf(val))
                arr[column * rows + row] = val;
    }

    void AddVal(int row, int column, double val)
    {
        if(arr[column * rows + row] + val < 0)
            return;
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

    void Normalize();
    void Brush(int centreX, int centreY, int radius, double val);
};


class Erosion
{
public:
    struct Config
    {
        int maxDropletLifetime = 120;
        double inertia = 0.4f;
        double erodeSpeed = 0.9f;
        double depositSpeed = 1.f;
        double sedimentCapacityFactor = 2.f;
        double evaporateSpeed = 0.05f;
        double minSlope = 0.0001f;
        double gravity = 10.f;
        double minSedimentCapacity = 0.01f;
        double stepFac = 1.f;
    };
    Grid *grid = nullptr;
    Config *config = nullptr;

private:
    int width, height;
    double initialWaterVolume = 1;
    double initialSpeed = 1;
    std::default_random_engine re;

private:
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

	double calcCarryCapacity(double heightDiff, double minSlope, double velocity, double water, double capacity);
    double calcHeightDiff(double HeightOld, double HeightNew);
    double calcVelocity(double velocity, double HeightDiff, double gravity);
    double calcWater(double water, double evaporation);
    gradient calcGradient(double posX, double posY);
    dropletMove CalcDropletMove(gradient g);
    dropletPoint CalcNewDropletPoint(double cx, double cy);
    void PlaceSediment(double amountToDeposit, int gridx, int gridy, double cellOffsetX, double cellOffsetY, double min,
                       double max, double current);

public:
    void DrawPath(sf::RenderWindow& window, double mx, double my, int cellsize);
    void Erode(double mx, double my);
    void Erode(int iter);
    Erosion(int width, int height, int iter, Config config);
};

