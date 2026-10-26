#include "demos.h"

#include <iostream>
#include <exception>

int main() {
    try {
        std::cout << "==============================================" << std::endl;
        std::cout << " АППРОКСИМАЦИЯ ДАННЫХ (Лекция 5)" << std::endl;
        std::cout << "==============================================" << std::endl << std::endl;

        demoExample8();
        demoExample9();
        demoExponential();
        demoExample10();
        demoExample11();

        std::cout << "==============================================" << std::endl;
        std::cout << " Все демонстрационные примеры выполнены" << std::endl;
        std::cout << "==============================================" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}