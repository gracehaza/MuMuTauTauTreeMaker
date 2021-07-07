import FWCore.ParameterSet.Config as cms

def customizeJets(process,coll,srcLabel='jets',postfix='',**kwargs):
    '''Customize jets'''
    isMC = kwargs.pop('isMC',False)
    reHLT = kwargs.pop('reHLT',False)
    jSrc = coll[srcLabel]
    rhoSrc = coll['rho']
    pvSrc = coll['vertices']

    # customization path
    pathName = 'jetCustomization{0}'.format(postfix)
    setattr(process,pathName,cms.Path())
    path = getattr(process,pathName)

    #################################
    ### add updated pileup jet id ###
    #################################
    # TODO: why is this here?
    #process.load("RecoJets.JetProducers.PileupJetID_cfi")
    #module = process.pileupJetId.clone(
    #    jets=cms.InputTag(jSrc),
    #    inputIsCorrected=True,
    #    applyJec=True,
    #    vertexes=cms.InputTag(pvSrc),
    #)
    #modName = 'pileupJetIdUpdated{0}'.format(postfix)
    #setattr(process,modName,module)

    #path *= getattr(process,modName)

    ######################
    ### recorrect jets ###
    ######################
    # TODO: reenable when we have a new recipe
    #from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

    #jetCorr = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')
    #if isMC:
    #    jetCorr = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')
    #updateJetCollection(
    #    process,
    #    jetSource = cms.InputTag(jSrc),
    #    jetCorrections = jetCorr,
    #    postfix=postfix,
    #)
    #modName = 'updatedPatJets{0}'.format(postfix)
    #getattr(process,modName).userData.userFloats.src += ['pileupJetIdUpdated{0}:fullDiscriminant'.format(postfix)]
    #jSrc = modName

    #################
    ### embed ids ###
    #################
    module = cms.EDProducer('DeepDiTauProducer',
        src = cms.InputTag(jSrc),
        DeepDiTauConfiguration = cms.PSet(
            memmapped = cms.bool(False),
            graphDefinitions = cms.VPSet(
                cms.PSet(
                    name = cms.string('ditau2017v1'),
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/03Dec2020_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/03Dec2020_mean_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('ditau2017MDv1'),
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/17Feb2021_relu_tanh_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/17Feb2021_relu_tanh_means_sigmas.txt'),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_nolepton_charm')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_charm_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_charm_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_nolepton_charm_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_charm_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_charm_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_nolepton_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_boosted_nolepton')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_boosted_nolepton_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_nolepton_charm_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_charm_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_charm_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_nolepton_charm')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_charm_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_charm_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_nolepton_massdeco')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_massdeco_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_massdeco_means_sigmas.txt'),
                ),
                cms.PSet(
                    name = cms.string('DeepDiTau_nolepton')
                    path = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_constantgraph.pb'),
                    means = cms.FileInPath('MuMuTauTauTreeMaker/MuTauTreelizer/data/June2021_nolepton_means_sigmas.txt'),
                ),
            ),
        ),
    )
    modName = 'deepDiTau{0}'.format(postfix)
    setattr(process,modName,module)
    path *= getattr(process,modName)

    module = cms.EDProducer(
        "JetIdEmbedder",
        src = cms.InputTag(jSrc),
        discriminator = cms.string('pileupJetId:fullDiscriminant'),
        ditau2017v1 = cms.InputTag("deepDiTau"+postfix,"ditau2017v1"),
        ditau2017MDv1 = cms.InputTag("deepDiTau"+postfix,"ditau2017MDv1"),
        DeepDiTau_boosted = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_boosted"),
        DeepDiTau_boosted_massdeco = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_boosted_massdeco"),
        DeepDiTau_boosted_nolepton_charm = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_boosted_nolepton_charm"),
        DeepDiTau_boosted_nolepton_charm_massdeco = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_boosted_nolepton_charm_massdeco"),
        DeepDiTau_boosted_nolepton_massdeco = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_boosted_nolepton_massdeco"),
        DeepDiTau_boosted_nolepton = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_boosted_nolepton"),
        DeepDiTau_massdeco = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_massdeco"),
        DeepDiTauValue = cms.InputTag("deepDiTau"+postfix,"DeepDiTauValue"),
        DeepDiTau_nolepton_charm_massdecoValue = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_nolepton_charm_massdeco"),
        DeepDiTau_nolepton_charmValue = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_nolepton_charm"),
        DeepDiTau_nolepton_massdecoValue = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_nolepton_massdeco"),
        DeepDiTau_noleptonValue = cms.InputTag("deepDiTau"+postfix,"DeepDiTau_nolepton"),
    )
    modName = 'jID{0}'.format(postfix)
    setattr(process,modName,module)
    jSrc = modName

    path *= getattr(process,modName)

    ###################
    ### embed truth ###
    ###################
    module = cms.EDProducer(
        "JetMCTruthEmbedder",
        src = cms.InputTag(jSrc),
        genSrc = cms.InputTag('prunedGenParticles'),
        packedGenSrc = cms.InputTag('packedGenParticles'),
    )
    modName = 'jTruth{0}'.format(postfix)
    setattr(process,modName,module)
    jSrc = modName

    path *= getattr(process,modName)

    #################
    ### embed rho ###
    #################
    module = cms.EDProducer(
        "JetRhoEmbedder",
        src = cms.InputTag(jSrc),
        rhoSrc = cms.InputTag(rhoSrc),
        label = cms.string("rho"),
    )
    modName = 'jRho{0}'.format(postfix)
    setattr(process,modName,module)
    jSrc = modName

    path *= getattr(process,modName)

    ##########################
    ### embed jet gen jets ###
    ##########################
    if isMC:
        module = cms.EDProducer(
            "JetGenJetEmbedder",
            src = cms.InputTag(jSrc),
            genJets = cms.InputTag("slimmedGenJets"),
            excludeLeptons = cms.bool(False),
            deltaR = cms.double(0.5),
        )
        modName = 'jGenJetMatching{0}'.format(postfix)
        setattr(process,modName,module)
        jSrc = modName

        path *= getattr(process,modName)

    ##############################
    ### embed trigger matching ###
    ##############################
    labels = []
    paths = []
    from triggers import triggerMap
    for trigger in triggerMap:
        if 'jet' in triggerMap[trigger]['objects']:
            labels += ['matches_{0}'.format(trigger)]
            paths += [triggerMap[trigger]['path']]
    module = cms.EDProducer(
        "JetHLTMatchEmbedder",
        src = cms.InputTag(jSrc),
        #triggerResults = cms.InputTag('TriggerResults', '', 'HLT'),
        triggerResults = cms.InputTag('TriggerResults', '', 'HLT2') if reHLT else cms.InputTag('TriggerResults', '', 'HLT'),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        deltaR = cms.double(0.4),
        labels = cms.vstring(*labels),
        paths = cms.vstring(*paths),
    )
    modName = 'jTrig{0}'.format(postfix)
    setattr(process,modName,module)
    jSrc = modName

    path *= getattr(process,modName)

    # add to schedule
    process.schedule.append(path)

    coll[srcLabel] = jSrc

    return coll
