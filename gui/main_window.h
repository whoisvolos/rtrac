#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_main_window.h"

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0, Qt::WindowFlags flags = 0);
    ~MainWindow();

private:
    Ui::MainWindow ui;

private slots:
    void openModelClicked();
};