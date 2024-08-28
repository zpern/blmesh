#ifndef CREATESOL_H
#define CREATESOL_H



namespace AABBSER
{
class bkMeshImpl;

/**
 * @brief 背景网格类，提供尺寸场查询功能
 */
class bkMesh
{
public:
    /**
     * @brief createFromBoundData   通过给定边界面网格输入，生成背景网格
     * @param bndPtNum      边界点数量
     * @param bndPts        边界点坐标 DIM = 3 * bndPtNum
     * @param bndPtSizes    边界点尺寸值 DIM = 1 * bndPtNum
     * @param bndFctNum     边界三角片数量
     * @param bndFcts       边界三角片索引值 DIM = 3 * bndFctNum
     * @param beta          尺寸随距离增长速率，需大于1
     * @return  成功则返回指向背景网格类实例的指针，请使用delete析构，失败则返回nullptr
     */
    static bkMesh* createFromBoundData(int bndPtNum, double bndPts[], double bndPtSizes[], int bndFctNum, int bndFcts[], double beta=1.2);

    /**
     * @brief getSize 输入坐标，返回该点在边界上投影点的尺寸值
     */
    double getSize(double x,double y,double z);

    /**
     * @brief getSizeSpatial 查询尺寸值
     * @return 返回该点在空间中的的尺寸值，失败则返回负数
     */
    double getSizeSpatial(double x,double y,double z);

    ~bkMesh();
private:
    bkMeshImpl* impl;
};

/**
 * @brief 使尺寸场函数兼容c风格函数指针的接口类
 * @details 使用方法:
 *          SizingFuncHelper<1>::setSingleton(new bkMesh(...))  //赋予背景网格对象
 *          SizingFuncHelper<1>::getSizingFunction()  //获取尺寸场函数指针
 */
template <int I>
class SizingFuncHelper{
public:
    typedef double (*SizingFunc_t)(double x,double y,double z);

    static bkMesh* getSingleton(){
        return singleton_;
    }

    static void setSingleton(bkMesh* singleton){
        singleton_ = singleton;
    }

    static SizingFunc_t getSizingFunction(){
        return &SizingFuncHelper::staticSizingFunction;
    }

private:
    static double staticSizingFunction(double x,double y,double z){
        return singleton_->getSizeSpatial(x,y,z);
    }

    static bkMesh* singleton_;
};

template <int I>
bkMesh* SizingFuncHelper<I>::singleton_ = nullptr;
}

#endif // CREATESOL_H
