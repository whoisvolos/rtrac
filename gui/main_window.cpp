#include <QFileDialog>
#include "gui/main_window.h"
#include <string>

MainWindow::~MainWindow() {

}

MainWindow::MainWindow(QWidget *parent, Qt::WindowFlags flags) : QMainWindow(parent, flags) {
    ui.setupUi(this);

    ui.centralwidget->setContentsMargins(6, 6, 6, 6);

    // open model
    connect(ui.actionOpen_obj, SIGNAL(triggered()), this, SLOT(openModelClicked()));
}

void MainWindow::openModelClicked() {

    QFileDialog* fileDialog = new QFileDialog(this);

    fileDialog->setFileMode(QFileDialog::ExistingFile);
    fileDialog->setNameFilter(QString::fromLocal8Bit("Model (*.obj)"));

    QStringList fileNames;

    if (fileDialog->exec()) {
        fileNames = fileDialog->selectedFiles();
        ui.statusbar->showMessage(fileNames.first(), 0);
    }

    /*
    if (fileNames.size() != 0) { // model loading
        string filePath = string(fileNames.at(0).toLocal8Bit().constData());

        // parsing model
        if (objLoader.LoadModel(filePath.c_str()) != -1) {
            string fileName = getFileNameFromPath(filePath);
            newModelLoadedInterfaceChanges(fileName);
        }
    }
    */

    delete fileDialog;
}
