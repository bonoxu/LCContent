/**
 *  @file   LCContent/src/LCMonitoring/EfficiencyMonitoringAlgorithm.cc
 * 
 *  @brief  Implementation of the efficiency monitoring algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCHelpers/SortingHelper.h"

#include "LCMonitoring/EfficiencyMonitoringAlgorithm.h"

using namespace pandora;

namespace lc_content
{

EfficiencyMonitoringAlgorithm::EfficiencyMonitoringAlgorithm() :
    m_mcThresholdEnergy(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EfficiencyMonitoringAlgorithm::~EfficiencyMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "sel", m_monitoringFileName, "RECREATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EfficiencyMonitoringAlgorithm::Run()
{

    // Extract the mc particle - this algorithm is intended to work only with single particle samples
    const MCParticleList *pMCParticleList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    if (1 != pMCParticleList->size())
    {
        std::cout << "EfficiencyMonitoring - Algorithm works only with single particle samples, nParticles " << pMCParticleList->size() << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    const MCParticle *const pMCParticle(*(pMCParticleList->begin()));

    // Extract the mc particle properties
    const float mcEnergy(pMCParticle->GetEnergy());
    const int mcPDGCode(pMCParticle->GetParticleId());
    const int mcPfoTarget(pMCParticle->IsPfoTarget());
    
    if (mcEnergy < m_mcThresholdEnergy)
    {
        std::cout << "EfficiencyMonitoring - MC particle energy below threshold " << mcEnergy << "( < " << m_mcThresholdEnergy << ")" <<std::endl;
        return STATUS_CODE_SUCCESS;
    }

    float radius(0.f), phi(0.f), theta(0.f);
    const CartesianVector &mcPosition(pMCParticle->GetEndpoint());
    mcPosition.GetSphericalCoordinates(radius, phi, theta);

    // Extract the most energetic pfo
    const PfoList *pPfoList= NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

    PfoVector pfoVector(pPfoList->begin(), pPfoList->end());
    std::sort(pfoVector.begin(), pfoVector.end(), SortingHelper::SortPfosByEnergy);
    const ParticleFlowObject *const pMostEnergeticPfo((!pfoVector.empty()) ? *(pfoVector.begin()) : NULL);

    // Extract the pfo properties
    float recoEnergy((pMostEnergeticPfo != NULL) ? pMostEnergeticPfo->GetEnergy() : 0.f);
    int recoPDGCode((pMostEnergeticPfo != NULL) ? pMostEnergeticPfo->GetParticleId() : 0);

    // Look for specific case of photon conversion to e+e-
    int isPhotonConversion(0);

    if ((mcPDGCode == PHOTON) && (std::abs(recoPDGCode) == E_MINUS) && (pfoVector.size() == 2))
    {
        const ParticleFlowObject *const pPfo1(pfoVector.at(0));
        const ParticleFlowObject *const pPfo2(pfoVector.at(1));

        if ((pPfo1->GetParticleId() == -(pPfo2->GetParticleId())) && (std::abs(pPfo1->GetParticleId()) == E_MINUS))
        {
            recoPDGCode = PHOTON;
            recoEnergy = pPfo1->GetEnergy() + pPfo2->GetEnergy();
            isPhotonConversion = 1;
        }
    }
    
    // test efficiency; use definition from CLIC particle ID efficiency
    int nMatchPFO(0), nPDGMatch(0), nEnergyMatch(0), nDirectionMatch(0);
    for (PfoVector::const_iterator iter = pfoVector.begin(), iterEnd = pfoVector.end(); iter != iterEnd; ++iter)
    {
        bool matchPFO(false), energyMatched(false), withinCone(false);
        if (pMCParticle && pMostEnergeticPfo)
        {
            if (0 != pMostEnergeticPfo->GetCharge())
            {
                const float mcPt2(pMCParticle->GetMomentum().GetX() * pMCParticle->GetMomentum().GetX() + pMCParticle->GetMomentum().GetY() * pMCParticle->GetMomentum().GetY());
                const float mcPt(mcPt2 > std::numeric_limits<float>::epsilon() ? std::sqrt(mcPt2) : std::numeric_limits<float>::epsilon());
                const float recoPt2(pMostEnergeticPfo->GetMomentum().GetX() * pMostEnergeticPfo->GetMomentum().GetX() + pMostEnergeticPfo->GetMomentum().GetY() * pMostEnergeticPfo->GetMomentum().GetY());
                const float recoPt(recoPt2 > std::numeric_limits<float>::epsilon() ? std::sqrt(recoPt2) : std::numeric_limits<float>::epsilon());
                energyMatched = std::fabs(mcPt - recoPt) < (0.05 * mcPt2);
                
                withinCone = std::fabs(pMCParticle->GetMomentum().GetOpeningAngle(pMostEnergeticPfo->GetMomentum())) < 0.0174533;
            }
            else
            {
                const float mcEroot(mcEnergy >  std::numeric_limits<float>::epsilon() ? std::sqrt(mcEnergy) : std::numeric_limits<float>::epsilon());
                energyMatched = std::fabs(mcEnergy - recoEnergy) < (4 * mcEroot + 0.5);
                withinCone = std::fabs(pMCParticle->GetMomentum().GetOpeningAngle(pMostEnergeticPfo->GetMomentum())) < 0.174533;
            }
            if (energyMatched && withinCone && (std::fabs(mcPDGCode) == std::fabs(recoPDGCode)))
                matchPFO = true;
        }
        if (energyMatched)
            nEnergyMatch++;
        if (withinCone)
            nDirectionMatch++;
        if (std::fabs(mcPDGCode) == std::fabs(recoPDGCode))
            nPDGMatch++;
        if (matchPFO)
            nMatchPFO++;
    }


    const float purity( (pfoVector.size() > 0) ? static_cast<double>(nMatchPFO) / static_cast<double>(pfoVector.size()) : 0.f);
    // Fill tree with information for this single particle event
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "mcPDGCode", mcPDGCode));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "recoPDGCode", recoPDGCode));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "mcEnergy", mcEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "mcPfoTarget", mcPfoTarget));
    
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "recoEnergy", recoEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "radius", radius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "phi", phi));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "theta", theta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "halfTheta", theta < std::atan(1)*2 ? theta : std::atan(1)*4 - theta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "isPhotonConversion", isPhotonConversion));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "PFOMatched", nMatchPFO));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "energyMatched", nEnergyMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "angleMatched", nDirectionMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "pdgMatched", nPDGMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "sel", "purity", purity));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "sel"));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EfficiencyMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MonitoringFileName", m_monitoringFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCThresholdEnergy", m_mcThresholdEnergy));

    return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
