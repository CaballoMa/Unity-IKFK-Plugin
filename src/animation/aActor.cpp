#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	m_BehaviorController = new BehaviorController();
	m_BehaviorController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_BehaviorController;
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

BehaviorController* AActor::getBehaviorController()
{
	return m_BehaviorController;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	if (!m_pSkeleton->getRootNode()) { return; }

	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	// 2.	Set the y component of the guide position to 0
	// 3.	Set the global rotation of the guide joint towards the guideTarget
	// Hint: Return of getGlobal***() here is in target space but not in world/global space
	vec3 oriPos = m_Guide.getGlobalTranslation();
	vec3 root = m_Guide.getLocal2Global().RotTrans(m_pSkeleton->getRootNode()->getGlobalTranslation());
	root[1] = 0;
	m_Guide.setGlobalTranslation(root);


	vec3 dir = guideTargetPos - m_Guide.getGlobalTranslation();
	dir[1] = 0;
	dir = dir.Normalize();

	vec3 oriDir = (m_Guide.getGlobalTranslation() - oriPos);
	oriDir[1] = 0;
	oriDir = oriDir.Normalize();
	mat3 rot = mat3::FromToRotation(oriDir, dir);
	m_Guide.setGlobalRotation(m_Guide.getGlobalRotation() * rot);
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { return; }
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space
	// Hint: Return of getGlobal***() here is in target space but not in world/global space

	// 1.	Update the local translation of the root based on the left height and the right height
	AJoint* root = m_pSkeleton->getRootNode();
	root->setLocalTranslation(root->getLocalTranslation() + vec3(0, (leftHeight + rightHeight) / 2, 0));
	m_pSkeleton->update();
	// 2.	Update the character with Limb-based IK 

	// Rotate Foot
	if (rotateLeft)
	{
		ATarget leftTarget;
		vec3 leftPos = leftFoot->getGlobalTranslation();
		leftTarget.setGlobalTranslation(vec3(leftPos[0], leftHeight, leftPos[2]));
		m_IKController->IKSolver_Limb(leftFoot->getID(), leftTarget);

		vec3 up = leftFoot->getGlobalRotation().Inverse() * leftNormal.Normalize();
		vec3 fwd = axisX.Cross(up).Normalize();
		vec3 bt = up.Cross(fwd).Normalize();
		mat3 rot = mat3::FromLocalAxis(bt, up, fwd);
		leftFoot->setLocalRotation(rot);
	}
	if (rotateRight)
	{
		ATarget rightTarget;
		vec3 rightPos = rightFoot->getGlobalTranslation();
		rightTarget.setGlobalTranslation(vec3(rightPos[0], rightHeight, rightPos[2]));
		m_IKController->IKSolver_Limb(rightFoot->getID(), rightTarget);

		vec3 up = rightFoot->getGlobalRotation().Inverse() * rightNormal.Normalize();
		vec3 fwd = axisX.Cross(up).Normalize();
		vec3 bt = up.Cross(fwd).Normalize();
		mat3 rot = mat3::FromLocalAxis(bt, up, fwd);
		rightFoot->setLocalRotation(rot);
	}
	m_pSkeleton->update();
}
