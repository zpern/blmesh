#include <spdlog/spdlog.h> 
 #include "../include/distancecalculation.h"
bool DistanceCalculator::isexist = false;
DistanceCalculator::DistanceCalculator()
{
	blm_ = (new BLMesh(cf.layer_num, 0, cf.isotropic_stop));
}

DistanceCalculator::~DistanceCalculator()
{
	isexist = false;
}
struct PoolGuard {
    bool prevFront, prevNode;
    PoolGuard(bool enable_pool) {
#ifdef USE_MEMORY_POOL
        prevFront = BLFront::PoolEnabled();
        prevNode  = BLNode::PoolEnabled();
        BLFront::SetPoolEnabled(enable_pool);
        BLNode::SetPoolEnabled(enable_pool);
#endif
    }
    ~PoolGuard() {
#ifdef USE_MEMORY_POOL
        BLFront::SetPoolEnabled(prevFront);
        BLNode::SetPoolEnabled(prevNode);
#endif
    }
};

void DistanceCalculator::ReadInput(INPUTFORMAT input,ConfigArgc cf)
{
	if (!isexist) {
		PoolGuard guard(false); 
		isexist = true;
        

		blm_->InitBLMesh();
		blm_->cf = cf;
        blm_->CalZeroNorm();
		blm_->cf.iscompresslen = false;
        blm_->generate_pyramid = true;
        blm_->SetBoundary(input);

		PreGenerate();
		std::cout.setf(ios::left, ios::adjustfield);
		std::cout.fill(' ');
		std::cout.setf(ios::left, ios::adjustfield);
		std::cout.fill(' ');
		//spdlog::info(std::setiosflags(ios::left) << std::setfill(' ') );
		spdlog::info("             Mesh Information    " );
		spdlog::info("+================================================+" );
		spdlog::info("|             Item                |number of cell|" );
		spdlog::info("-------------------------------------------------" );
		spdlog::info("|            Pyramid              |{}", blm_->m_nPyramid );
		spdlog::info("|            Tetrahedron          |{}", blm_->m_nTet );
		spdlog::info("|         Triangular prism        |{}", blm_->m_nPrism );
		spdlog::info("|            Triangle             |{}", blm_->m_nTri );
		spdlog::info("|          Quadrilateral          |{}", blm_->m_nQuad );

		Calculate();
	//	blm_->RemvNodElm();
	//	blm_->FreeMemoryInFrontAndNode();

		delete blm_;
	}
}

std::vector<double> DistanceCalculator::getHeightProjection()
{

	return recommand_array_;
}

void DistanceCalculator::PreGenerate()
{
	try {
		blm_->GenerateBLMesh();
	}
	catch (string error) {
		spdlog::info(error );
	}
}

void DistanceCalculator::Calculate()
{
	recommand_array_ = blm_->recommand_length_calculation();

}
